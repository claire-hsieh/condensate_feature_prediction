import os
import random
import pickle
import psutil
from pynvml import nvmlInit, nvmlDeviceGetHandleByIndex, nvmlDeviceGetTemperature, NVML_TEMPERATURE_GPU
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
from scipy.stats import pearsonr, spearmanr
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from sklearn.model_selection import KFold, LeaveOneGroupOut, GroupKFold, train_test_split
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay, accuracy_score, f1_score
from collections import Counter

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.distributions import Categorical
from torch.utils.data import DataLoader, Dataset, Subset, TensorDataset
import subprocess

####################
### Load in Data ###
####################

small_pl_data = pd.read_csv("/Users/clairehsieh/OneDrive/Documents/UCLA/rotations/kalli_kappel/data/small_pool_data_with_seq_info_202507.csv", 
                             index_col=False)
list(small_pl_data.columns)
small_pl_data.shape

# large pool data
large_pl_data = pd.read_csv("/Users/clairehsieh/OneDrive/Documents/UCLA/rotations/kalli_kappel/data/large_pool_data_with_seq_info_202507.csv")
large_pl_data.shape
list(large_pl_data.columns)

#####################
### Preprocessing ###
#####################

def read_mmseqs_into_distance(filename):
    df = pd.read_csv(filename, sep="\t", usecols=[0,1,2], names=["query", "target", "identity"])
    seq_ids = pd.unique(df[['query','target']].values.ravel())
    n = len(seq_ids)

    # Map sequence ID to index
    id_to_idx = {seq_id: i for i, seq_id in enumerate(seq_ids)}
    dist_matrix = np.ones((n,n))
    np.fill_diagonal(dist_matrix, 0)
    for _, row in df.iterrows():
        i = id_to_idx[row['query']]
        j = id_to_idx[row['target']]
        dist = 1.0 - row['identity']  # identity should be in [0,1]; divide by 100 if in %
        dist_matrix[i,j] = dist
        dist_matrix[j,i] = dist  # symmetric
    return dist_matrix

def hierarchical_cluster_mmseqs(filename, num_clusters=10, threshold=0.3, method='threshold'):
    if method not in ["threshold", "num_clusters"]: raise ValueError("method should be 'threshold' or 'num_clusters'")
    df = pd.read_csv(filename, sep="\t", usecols=[0,1,2], names=["query", "target", "identity"])
    seq_ids = pd.unique(df[['query','target']].values.ravel())
    dist_matrix = read_mmseqs_into_distance(filename)

    Z = linkage(dist_matrix, method='average')  # or 'single', 'complete', 'ward'
    dendrogram(Z, labels=seq_ids)
    plt.show()

    # Cut tree to get cluster labels    
    if method == "threshold": clusters = fcluster(Z, threshold, criterion='distance')
    elif method == "num_clusters": clusters = fcluster(Z, t=num_clusters, criterion='maxclust')
    return clusters

# one hot encode protein seq
def one_hot_encode(seq):
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    aa_to_index = {aa: i for i, aa in enumerate(amino_acids)}

    # shape: (length_of_seq, 20)
    encoding = np.zeros((len(seq), len(amino_acids)), dtype=np.int8)
    for i, aa in enumerate(seq):
        if aa in aa_to_index:  # skip if non-standard residue
            encoding[i, aa_to_index[aa]] = 1
    return encoding

# create kmer table
def kmer_table(seqs, k=3):
    """
    Convert a list of sequences into a k-mer count matrix.    
    Args:
        seqs (list of str): List of sequences (all same alphabet, e.g. DNA/protein).
        k (int): k-mer size.        
    Returns:
        np.ndarray: Matrix of shape (num_seqs, num_unique_kmers).
        list: List of unique kmers in column order.
    """

    # collect all kmers across all seqs
    all_kmers = []
    for seq in seqs:
        for i in range(len(seq) - k + 1):
            all_kmers.append(seq[i:i+k])
    unique_kmers = sorted(set(all_kmers))  # column order

    # mapping for column index
    kmer_index = {kmer: i for i, kmer in enumerate(unique_kmers)}

    # build matrix
    matrix = np.zeros((len(seqs), len(unique_kmers)), dtype=int)
    for row, seq in enumerate(seqs):
        counts = Counter(seq[i:i+k] for i in range(len(seq) - k + 1))
        for kmer, c in counts.items():
            col = kmer_index[kmer]
            matrix[row, col] = c

    return matrix, unique_kmers


################ 
### Modeling ###
################

class Config:
    def __init__(self,
                 output_dir,
                 seed=0,
                 num_folds=5,
                 max_epochs=100,
                 lr=0.001,
                 batch_size=10,
                 input_dim=100,
                 hidden_dim=512,
                 output_dim=1,
                 num_epochs=500,
                 num_layers=3,
                 verbose=False,
                 layernorm=False,
                 dropout=False,
                 dropout_rate=0.05, 
                 early_stopping=True,
                 patience = 30,
                 delta = 0.05, 
                 train_split = "kfold", 
                 criterion = "regression", 
                 save_model = False):       
        self.layernorm = layernorm
        self.dropout = dropout
        self.dropout_rate = dropout_rate
        if torch.cuda.is_available():            self.device = torch.device("cuda")
        elif torch.backends.mps.is_available():  self.device = torch.device("mps")
        else:                                    self.device = torch.device("cpu")
        print("Using device:", self.device)    
            
        # training parameters
        self.seed = seed
        self.num_folds = num_folds
        self.max_epochs = max_epochs
        self.lr = lr
        self.input_dim = input_dim
        self.batch_size = batch_size
        self.hidden_dim = hidden_dim
        self.output_dim = output_dim
        self.num_epochs = num_epochs
        self.num_layers = num_layers
        self.verbose = verbose
        self.early_stopping = early_stopping
        self.patience = patience    
        self.delta = delta
        self.output_dir = f"{output_dir}/seed_{seed}/lr_{lr}/hidden_dim_{hidden_dim}/num_layers_{num_layers}/batch_size_{batch_size}/"
        self.train_split = train_split
        self.criterion = criterion
        self.save_model = save_model
        if not os.path.exists(self.output_dir): os.makedirs(self.output_dir)
        if dropout: self.output_dir += f"dropout_{dropout_rate}/"
        if layernorm: self.output_dir += f"layernorm/"
        if early_stopping: self.output_dir += f"early_stopping_patience_{patience}_delta_{delta}/"

    def __str__(self):
        return '\n'.join(f'{k}: {v}' for k, v in vars(self).items())

    def __repr__(self):
        return self.__str__()

class MLP(nn.Module):
    def __init__(self, config):
        super().__init__()
        self.hidden_dim = config.hidden_dim
        self.num_layers = config.num_layers
        self.input_dim = config.input_dim
        
        self.mlp = nn.Sequential(
            nn.Linear(self.input_dim , self.hidden_dim),
            nn.LeakyReLU()
        ).to(config.device)
        
        # final layer always outputs 1 value per sample
        self.final = nn.Linear(self.hidden_dim, 1).to(config.device)
        
        self.layer_norm = nn.LayerNorm(self.hidden_dim).to(config.device)
        self.logZ = nn.Parameter(torch.ones(1, device=config.device))

    def forward(self, x):
        output = self.mlp(x)        # [batch_size, hidden_dim]
        output = self.final(output) # [batch_size, 1]
        return output.squeeze(-1)   # [batch_size]

    def _init_weights(self):
        for module in self.modules():
            if isinstance(module, nn.Linear):
                nn.init.xavier_uniform_(module.weight)
                if module.bias is not None:
                    nn.init.zeros_(module.bias)


class CondensateDataset(Dataset):
    def __init__(self, x, y):
        # x: 2D array-like [n_samples, n_features]
        # y: 1D array-like [n_samples]
        self.x = torch.tensor(x, dtype=torch.float32)
        self.y = torch.tensor(y, dtype=torch.float32)

    def __len__(self):
        return len(self.x)

    def __getitem__(self, idx):
        return self.x[idx], self.y[idx]

#### UTILS ##### 
class EarlyStopping:
    def __init__(self, patience=15, delta=0):
        self.patience = patience
        self.delta = delta
        self.best_score = None
        self.early_stop = False
        self.counter = 0
        self.best_model_state = None
    def __call__(self, val_loss, model):
        score = -val_loss
        if self.best_score is None:
            self.best_score = score
            self.best_model_state = model.state_dict()
        elif score < self.best_score + self.delta:
            self.counter += 1
            if self.counter >= self.patience:
                self.early_stop = True
        else:
            self.best_score = score
            self.best_model_state = model.state_dict()
            self.counter = 0
    def load_best_model(self, model):
        model.load_state_dict(self.best_model_state)

def log_memory_usage():
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    print(f"Memory Usage: {mem_info.rss / (1024 ** 2):.2f} MB")

def mean(ls):
    return sum(ls) / len(ls)

def compute_spearman_correlation(actuals, predictions):
    try:
        if np.all(actuals == actuals[0]) or np.all(predictions == predictions[0]):
            print("Warning: Constant array detected - correlation undefined")
            return 0          
        correlation, _ = spearmanr(actuals, predictions)
        return correlation
    except Exception as e:
        print(f"Error computing correlation: {e}")
        return None

def compute_pearson_correlation(actuals, predictions):
    if isinstance(actuals, torch.Tensor):
        actuals = actuals.cpu().numpy()
    if isinstance(predictions, torch.Tensor):
        predictions = predictions.cpu().numpy()
    correlation, _ = pearsonr(actuals, predictions)
    return correlation

def plot_epoch_losses(epoch_train_losses, epoch_val_losses, filename, output_dir):
    print(f"Train Losses: {epoch_train_losses}")
    print(f"Val Losses: {epoch_val_losses}")
    plt.figure()
    plt.plot(epoch_train_losses, label='Training Loss')
    plt.plot(epoch_val_losses, label='Validation Loss')
    plt.title(f'Training and Validation Loss | Fold {filename}')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.legend()
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    plt.savefig(output_dir + f"Fold{filename}_loss_plot.jpg")
    plt.close()

def plot_fold_scatter(val_actual, val_pred, train_actual, train_pred, correlations, fold, output_dir):   
    # Plot scatterplot for validation data
    plt.figure(figsize=(8, 6))
    plt.scatter(val_actual, val_pred, alpha=0.6, label="Validation Data", color="orange")
    plt.plot([min(val_actual), max(val_actual)], [min(val_actual), max(val_actual)], 'r--', label="Ideal Fit")
    plt.xlabel("Actual Values")
    plt.ylabel("Predicted Values")
    plt.title("Validation Data: Actual vs Predicted")
    plt.legend()
    # Add correlations as text on the validation plot
    plt.text(
        0.95, 0.05,
        f"Spearman: {correlations['val_spearman']:.2f}\nPearson: {correlations['val_pearson']:.2f}",
        transform=plt.gca().transAxes,
        fontsize=10,
        verticalalignment='bottom',
        horizontalalignment='right',
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.5)
    )
    plt.savefig(os.path.join(output_dir, f"validation_fold_{fold}.png"))
    plt.close()

    # Plot scatterplot for training data
    plt.figure(figsize=(8, 6))
    plt.scatter(train_actual, train_pred, alpha=0.6, label="Training Data", color="blue")
    plt.plot([min(train_actual), max(train_actual)], [min(train_actual), max(train_actual)], 'r--', label="Ideal Fit")
    plt.xlabel("Actual Values")
    plt.ylabel("Predicted Values")
    plt.title("Training Data: Actual vs Predicted")
    plt.legend()
     # Add correlations as text on the training plot
    plt.text(
        0.95, 0.05,
        f"Spearman: {correlations['train_spearman']:.2f}\nPearson: {correlations['train_pearson']:.2f}",
        transform=plt.gca().transAxes,
        fontsize=10,
        verticalalignment='top',
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.5)
    )
    plt.savefig(os.path.join(output_dir, f"training_fold_{fold}.jpg"))
    plt.close()

def plot_all_folds_losses(losses, output_dir, filename="all_folds"):
    """
    losses: dict with keys 'train_losses' and 'val_losses'
            values are lists of lists, one list per fold
            e.g. losses['train_losses'][i] = list of training losses for fold i
    output_dir: directory to save plot
    filename: filename prefix for saved plot
    """
    train_folds = losses["train_losses"]
    val_folds = losses["val_losses"]

    n_folds = len(train_folds)
    fig, axes = plt.subplots(1, n_folds, figsize=(5 * n_folds, 4), sharey=True)

    if n_folds == 1:  # make axes iterable if only one subplot
        axes = [axes]

    for i, ax in enumerate(axes):
        ax.plot(train_folds[i], label="Training Loss")
        ax.plot(val_folds[i], label="Validation Loss")
        ax.set_title(f"Fold {i+1}")
        ax.set_xlabel("Epoch")
        ax.set_ylabel("Loss")
        ax.legend()

    plt.tight_layout()
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    plt.savefig(os.path.join(output_dir, f"{filename}_loss_plot.jpg"))
    plt.close()

def plot_all_folds_scatter(fold_outputs, all_correlations, output_dir, filename="all_folds_scatter"):
    """
    fold_outputs: dict with keys
        'fold_train_actual', 'fold_train_pred',
        'fold_val_actual', 'fold_val_pred'
        Each is a list of lists (one per fold).

    all_correlations: list of dicts (one per fold), each dict has keys
        'train_spearman', 'val_spearman', 'train_pearson', 'val_pearson'

    output_dir: directory to save the figure
    filename: filename prefix for saved plot
    """

    n_folds = len(fold_outputs["fold_train_actual"])
    fig, axes = plt.subplots(2, n_folds, figsize=(5 * n_folds, 10), sharey=True)

    if n_folds == 1:  # axes will be 1D
        axes = axes[:, None]  # make it 2x1 array

    for i in range(n_folds):
        # Training subplot (top row)
        ax_train = axes[0, i]
        train_actual = fold_outputs["fold_train_actual"][i]
        train_pred = fold_outputs["fold_train_pred"][i]
        ax_train.scatter(train_actual, train_pred, alpha=0.6, color="blue", label="Training")
        ax_train.plot([min(train_actual), max(train_actual)],
                      [min(train_actual), max(train_actual)],
                      'r--', label="Ideal Fit")
        ax_train.set_xlabel("Actual Values")
        ax_train.set_ylabel("Predicted Values")
        ax_train.set_title(f"Fold {i+1} Training")
        ax_train.legend()
        # Correlation text
        corr_val = all_correlations[i]
        ax_train.text(
            0.95, 0.05,
            f"Spearman: {corr_val['train_spearman']:.2f}\nPearson: {corr_val['train_pearson']:.2f}",
            transform=ax_train.transAxes,
            fontsize=10,
            verticalalignment='bottom',
            horizontalalignment='right',
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.5)
        )

        # Validation subplot (bottom row)
        ax_val = axes[1, i]
        val_actual = fold_outputs["fold_val_actual"][i]
        val_pred = fold_outputs["fold_val_pred"][i]
        ax_val.scatter(val_actual, val_pred, alpha=0.6, color="orange", label="Validation")
        ax_val.plot([min(val_actual), max(val_actual)],
                    [min(val_actual), max(val_actual)],
                    'r--', label="Ideal Fit")
        ax_val.set_xlabel("Actual Values")
        ax_val.set_ylabel("Predicted Values")
        ax_val.set_title(f"Fold {i+1} Validation")
        ax_val.legend()
        # Correlation text
        ax_val.text(
            0.95, 0.05,
            f"Spearman: {corr_val['val_spearman']:.2f}\nPearson: {corr_val['val_pearson']:.2f}",
            transform=ax_val.transAxes,
            fontsize=10,
            verticalalignment='bottom',
            horizontalalignment='right',
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.5)
        )

    plt.tight_layout()
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    plt.savefig(os.path.join(output_dir, f"{filename}.png"))
    plt.close()

def get_accuracy_f1(actuals_val, preds_val, actuals_train, preds_train):
    # Compute accuracy and F1 score
    train_accuracy = accuracy_score(actuals_train, preds_train)
    val_accuracy = accuracy_score(actuals_val, preds_val)
    train_f1 = f1_score(actuals_train, preds_train, average='weighted')
    val_f1 = f1_score(actuals_val, preds_val, average='weighted')
    return train_accuracy, val_accuracy, train_f1, val_f1

def get_baseline_accuracy(actuals_val, actuals_train):
    # ensure integer labels (0/1)
    actuals_train = np.rint(actuals_train).astype(int)
    actuals_val   = np.rint(actuals_val).astype(int)

    majority_class_train = np.argmax(np.bincount(actuals_train))
    baseline_accuracy_train = np.mean(actuals_train == majority_class_train)

    majority_class_val = np.argmax(np.bincount(actuals_val))
    baseline_accuracy_val = np.mean(actuals_val == majority_class_val)

    return baseline_accuracy_train, baseline_accuracy_val

def plot_confusion_matrix(fold_outputs, config, name=None):
    n_folds = len(fold_outputs["fold_train_actual"])
    averages = {"train_baseline":[], "train": [], "val_baseline":[], "val": []}

    # 2 rows (train, val), n_folds columns
    fig, axes = plt.subplots(2, n_folds, figsize=(7 * n_folds, 12))

    if n_folds == 1:  # ensure consistent indexing
        axes = axes.reshape(2, 1)

    for fold in range(n_folds):
        actuals_val   = np.rint(fold_outputs["fold_val_actual"][fold]).astype(int)
        preds_val     = np.rint(fold_outputs["fold_val_pred"][fold]).astype(int)
        actuals_train = np.rint(fold_outputs["fold_train_actual"][fold]).astype(int)
        preds_train   = np.rint(fold_outputs["fold_train_pred"][fold]).astype(int)

        cm_train = confusion_matrix(actuals_train, preds_train)
        cm_val   = confusion_matrix(actuals_val, preds_val)

        baseline_accuracy_train, baseline_accuracy_val = get_baseline_accuracy(actuals_val, actuals_train)
        train_accuracy, val_accuracy, train_f1, val_f1 = get_accuracy_f1(
            actuals_val, preds_val, actuals_train, preds_train
        )

        averages["train_baseline"].append(baseline_accuracy_train)
        averages["val_baseline"].append(baseline_accuracy_val)
        averages["train"].append(train_accuracy)
        averages["val"].append(val_accuracy)

        # Train subplot (row 0)
        ax_train = axes[0, fold]
        sns.heatmap(cm_train, annot=True, fmt='d', cmap='Blues', annot_kws={"size": 12}, ax=ax_train)
        ax_train.set_title(
            f'Training (Fold {fold})\nAcc: {train_accuracy:.2f}, F1: {train_f1:.2f}, Baseline: {baseline_accuracy_train:.2f}',
            fontsize=14
        )
        ax_train.set_xlabel('Predicted', fontsize=12)
        ax_train.set_ylabel('Actual', fontsize=12)

        # Validation subplot (row 1)
        ax_val = axes[1, fold]
        sns.heatmap(cm_val, annot=True, fmt='d', cmap='Reds', annot_kws={"size": 12}, ax=ax_val)
        ax_val.set_title(
            f'Validation (Fold {fold})\nAcc: {val_accuracy:.2f}, F1: {val_f1:.2f}, Baseline: {baseline_accuracy_val:.2f}',
            fontsize=14
        )
        ax_val.set_xlabel('Predicted', fontsize=12)
        ax_val.set_ylabel('Actual', fontsize=12)

    # averages summary below plot
    avg_train_acc  = np.mean(averages["train"])
    avg_val_acc    = np.mean(averages["val"])
    avg_train_base = np.mean(averages["train_baseline"])
    avg_val_base   = np.mean(averages["val_baseline"])

    fig.text(0.5, 0.01, 
             f"Averages \n Train Acc: {avg_train_acc:.2f}, Val Acc: {avg_val_acc:.2f}\n "
             f"Train Baseline: {avg_train_base:.2f}, Val Baseline: {avg_val_base:.2f}",
             ha='center', fontsize=14)

    plt.tight_layout(rect=[0, 0.05, 1, 1])
    plt.savefig(os.path.join(config.output_dir, f'confusion_matrices_all_folds.png'),
                bbox_inches="tight")
    # plt.show()
    plt.close()

def get_predictions(model, dataloader, config):
    model.eval()
    all_actuals = []
    all_predictions = []
    with torch.no_grad():
        for batch in dataloader:
            x,y = batch
            output = model(x)
            if config.criterion == "classification": 
                output = (torch.sigmoid(torch.tensor(output)) >= 0.5).int() #logits --> probabilities --> class]
            all_actuals.extend(y.cpu().numpy())
            all_predictions.extend(output.cpu().numpy())
    all_actuals = np.array(all_actuals)
    all_predictions = np.array(all_predictions)
    return all_actuals, all_predictions

def plot_prediction_results(actual, predicted, filename = None):
    plt.scatter(actual, predicted)
    plt.xlabel("actual")
    plt.ylabel("predicted")
    plt.title("Actual vs Predicted")
    plt.show()
    if filename: plt.savefig(filename)

def plot_test_scatter(test_actual, test_pred, config):
    plt.figure(figsize=(6,6))
    plt.scatter(test_actual, test_pred, alpha=0.5)
    plt.plot([min(test_actual), max(test_actual)], [min(test_actual), max(test_actual)], 'r--')
    plt.xlabel("Actual")
    plt.ylabel("Predicted")
    plt.title("Test Set Scatter")
    plt.savefig(os.path.join(config.output_dir, "test_scatter.png"))
    plt.close()
    print(f"Test set scatterplot saved to {os.path.join(config.output_dir, 'test_scatter.png')}")

def plot_test_confusion_matrix(test_actual, test_pred, config):
    actuals = np.rint(np.array(test_actual)).astype(int)
    preds   = np.rint(np.array(test_pred)).astype(int)
    baseline_accuracy = max(np.mean(actuals == 0), np.mean(actuals == 1))  # majority class
    accuracy = accuracy_score(actuals, preds)
    cm = confusion_matrix(actuals, preds)
    plt.figure(figsize=(6,6))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Reds')
    plt.xlabel("Predicted")
    plt.ylabel("Actual")
    plt.title(f"Test Confusion Matrix\nAccuracy: {accuracy:.2f}, Baseline: {baseline_accuracy:.2f}")
    plt.savefig(os.path.join(config.output_dir, "test_confusion_matrix.png"))
    plt.close()
    print(f"Test set cm saved to {os.path.join(config.output_dir, 'test_confusion_matrix.png')}")

def train(X, y, config, groups=None):    
    dataset = CondensateDataset(X, y)

    output_dir = config.output_dir
    fold_outputs = {"fold_train_actual":[], 
                    "fold_train_pred":[],
                     "fold_val_actual":[], 
                     "fold_val_pred":[]}
    logs = {"fold_actuals":[],
            "fold_preds":[],
            "fold_train_actuals":[],
            "fold_train_preds":[],
            "fold_train_losses":[],
            "fold_val_losses":[],
            "train_num_genes":[],
            "all_train_labels":[],
            "val_num_genes":[],
            "all_val_labels":[],
            "f1_accuracy":{},
            "all_correlation":[]}
    if config.criterion == "regression":     loss_fn = nn.MSELoss()
    elif config.criterion == "classification": loss_fn = nn.BCEWithLogitsLoss()
    else: print("config.criterion should be 'regression' or 'classification'")

    training_logs = {"train_actual" : [],
                    "train_pred" : [],
                    "val_actual" : [],
                    "val_pred" : [],
                    "epoch_train_losses" : [],
                    "epoch_val_losses" : []}
    model_output_dir = output_dir + "/model/"
    if not os.path.exists(model_output_dir): os.makedirs(model_output_dir)

    losses = {"train_losses":[], "val_losses":[]}
    all_correlations = []

    #### Training loop ##### 

    #  Train test split   
    X_train, X_test, y_train, y_test, groups_train, groups_test = train_test_split(
        X.detach().cpu().numpy(), y.detach().cpu().numpy(), groups,
        test_size=0.15,
        random_state=config.seed,
        stratify=y.detach().cpu().numpy() if config.criterion == "classification" else None
    )

    X_train = torch.tensor(X_train, dtype=torch.float32, device=config.device)
    X_test  = torch.tensor(X_test, dtype=torch.float32, device=config.device)
    y_train = torch.tensor(y_train, dtype=torch.float32, device=config.device)
    y_test  = torch.tensor(y_test, dtype=torch.float32, device=config.device)

    train_dataset = CondensateDataset(X_train, y_train)
    test_dataset  = CondensateDataset(X_test, y_test)

    if config.train_split == "kfold":
        splitter = KFold(n_splits=config.num_folds, shuffle=True, random_state=config.seed)
        splits = list(splitter.split(range(len(train_dataset))))

    elif config.train_split == "group_kfold":
        splitter = GroupKFold(n_splits=config.num_folds)
        splits = list(splitter.split(X_train, y_train, groups=groups_train))

    else:
        splitter = LeaveOneGroupOut()
        splits = list(splitter.split(X_train, y_train, groups_train))
    

    for fold, (train_idx, val_idx) in enumerate(splits):
        print(f"Fold {fold+1}/{len(splits)}, train size={len(train_idx)}, val size={len(val_idx)}")
        model = MLP(config)
        train_dataloader = DataLoader(Subset(dataset, train_idx), batch_size=config.batch_size, shuffle=True, num_workers=0)
        val_dataloader = DataLoader(Subset(dataset, val_idx), batch_size=config.batch_size, shuffle=False, num_workers=0)
        epoch_train_losses = []
        epoch_val_losses = []
        val_labels = []
        # for fold in range(num_folds):
        filename = fold
        optimizer = optim.Adam(model.parameters(), lr=config.lr)
        if config.early_stopping:
            early_stopping =  EarlyStopping(patience=config.patience, delta=config.delta)
        print(f'Fold {fold+1}/{len(splits)}')
        log_memory_usage()

        for epoch in tqdm(range(config.max_epochs)):
            model.train()
            train_losses = []
            val_losses = []
            for batch in train_dataloader:
                x, y = batch
                x = x.to(config.device)
                y = y.to(config.device)
                optimizer.zero_grad()
                device = next(model.parameters()).device  
                output = model(x)
                if config.criterion == "regression":       loss = loss_fn(output.float(), y.float())
                elif config.criterion == "classification": loss = loss_fn(output.view(-1), y.float().view(-1))

                train_losses.append(loss.item())
                loss.backward()
                optimizer.step()
            model.eval()
            with torch.no_grad():
                for batch in val_dataloader:
                    x, y = batch
                    x = x.to(config.device)
                    y = y.to(config.device)
                    output = model(x)
                    loss = loss_fn(output.float(), y.float())                
                    val_losses.append(loss.item())
                epoch_train_losses.append(mean(train_losses))
                epoch_val_losses.append(mean(val_losses))
                if config.early_stopping:
                    early_stopping(mean(val_losses), model)
                    if early_stopping.early_stop:
                        print(f"Early stopping at epoch {epoch}")
                        early_stopping.load_best_model(model)
                        break
        
        logs["all_val_labels"].append(val_labels)
        logs["fold_train_losses"].append(epoch_train_losses)
        logs["fold_val_losses"].append(epoch_val_losses)
        train_actual, train_pred = get_predictions(model, train_dataloader, config)
        val_actual, val_pred = get_predictions(model, val_dataloader, config)
        logs["fold_actuals"].append(val_actual)
        logs["fold_preds"].append(val_pred)
        logs["fold_train_actuals"].append(train_actual)
        logs["fold_train_preds"].append(train_pred)
        if config.criterion == "regression": 
            correlations = {
                "train_spearman": compute_spearman_correlation(train_actual, train_pred) or 0,
                "val_spearman": compute_spearman_correlation(val_actual, val_pred) or 0,
                "train_pearson": compute_pearson_correlation(train_actual, train_pred) or 0,
                "val_pearson": compute_pearson_correlation(val_actual, val_pred) or 0
            }
            all_correlations.append(correlations)
            logs["all_correlation"].append(correlations)
        # plot_epoch_losses(epoch_train_losses, epoch_val_losses, filename, output_dir)
        # plot_fold_scatter(val_actual, val_pred, train_actual, train_pred, correlations, fold, output_dir)        

        losses["train_losses"].append(epoch_train_losses), losses["val_losses"].append(epoch_val_losses)
        fold_outputs["fold_train_actual"].append(train_actual), fold_outputs["fold_train_pred"].append(train_pred), fold_outputs["fold_val_actual"].append(val_actual), fold_outputs["fold_val_pred"].append(val_pred), 
        log_memory_usage()
        if config.save_model: torch.save(model.state_dict(), os.path.join(model_output_dir, f"model_fold_{fold}.pth"))
        # plot_correlations(logs["all_correlation"], output_dir)

    # Test set
    print("Evaluating test set")
    test_loader = DataLoader(test_dataset, batch_size=config.batch_size, shuffle=False)
    test_actual, test_pred = get_predictions(model, test_loader, config)
    if config.criterion == "regression":
        test_metrics = {
            "test_spearman": compute_spearman_correlation(test_actual, test_pred),
            "test_pearson": compute_pearson_correlation(test_actual, test_pred),
        }
        plot_all_folds_scatter(fold_outputs, all_correlations, config.output_dir, filename="all_folds_scatter")
        plot_test_scatter(test_actual, test_pred, config)
    elif config.criterion == "classification":
        test_metrics = {"test_confusion": confusion_matrix(test_actual, (torch.sigmoid(torch.tensor(test_pred)) > 0.5).int())}
        plot_confusion_matrix(fold_outputs, config, name=None)
        plot_test_confusion_matrix(test_actual, test_pred, config)

    logs["test_results"] = test_metrics

    # Save logs and fold_outputs as pkl
    plot_all_folds_losses(losses, output_dir, filename="all_folds")
    if config.criterion == "regression":     plot_all_folds_scatter(fold_outputs, all_correlations, output_dir, filename="all_folds_scatter")
    elif config.criterion == "classification": plot_confusion_matrix(fold_outputs, config, name=None)

    logs_file = os.path.join(output_dir, "logs.pkl")
    with open(logs_file, "wb") as f:
        pickle.dump(logs, f)
    print(f"Logs outputs saved to {logs_file}")
    fold_outputs_file = os.path.join(output_dir, "fold_outputs.pkl")
    with open(fold_outputs_file, "wb") as f:
        pickle.dump(fold_outputs, f)
    print(f"Fold outputs saved to {fold_outputs_file}")

    return model, logs, fold_outputs, losses

def train_with_predefined_splits(model_input, config):    
    output_dir = config.output_dir
    fold_outputs = {"fold_train_actual":[], 
                    "fold_train_pred":[],
                     "fold_val_actual":[], 
                     "fold_val_pred":[]}
    logs = {"fold_actuals":[],
            "fold_preds":[],
            "fold_train_actuals":[],
            "fold_train_preds":[],
            "fold_train_losses":[],
            "fold_val_losses":[],
            "train_num_genes":[],
            "all_train_labels":[],
            "val_num_genes":[],
            "all_val_labels":[],
            "f1_accuracy":{},
            "all_correlation":[]}
    if config.criterion == "regression":     loss_fn = nn.MSELoss()
    elif config.criterion == "classification": loss_fn = nn.BCEWithLogitsLoss()
    else: print("config.criterion should be 'regression' or 'classification'")

    training_logs = {"train_actual" : [],
                    "train_pred" : [],
                    "val_actual" : [],
                    "val_pred" : [],
                    "epoch_train_losses" : [],
                    "epoch_val_losses" : []}
    model_output_dir = output_dir + "/model/"
    if not os.path.exists(model_output_dir): os.makedirs(model_output_dir)

    losses = {"train_losses":[], "val_losses":[]}
    all_correlations = []

    #### Training loop ##### 

    #  Train test split   
    train_dataset = CondensateDataset(model_input["train"]["x"], model_input["train"]["y"])
    val_dataset   = CondensateDataset(model_input["val"]["x"], model_input["val"]["y"])
    test_dataset  = CondensateDataset(model_input["test"]["x"],model_input["test"]["y"])

    # K-fold on training dataset
    kf = KFold(n_splits=5, shuffle=True, random_state=42)

    for fold, (train_idx, val_idx) in enumerate(kf.split(range(len(train_dataset)))):
        print(f"Fold {fold+1}")
        train_dataloader = DataLoader(Subset(train_dataset, train_idx), batch_size=config.batch_size, shuffle=True)
        val_dataloader = DataLoader(Subset(train_dataset, val_idx), batch_size=config.batch_size, shuffle=False)


        model = MLP(config)
        epoch_train_losses = []
        epoch_val_losses = []
        val_labels = []
        # for fold in range(num_folds):
        filename = fold
        optimizer = optim.Adam(model.parameters(), lr=config.lr)
        if config.early_stopping:
            early_stopping =  EarlyStopping(patience=config.patience, delta=config.delta)
        log_memory_usage()

        for epoch in tqdm(range(config.max_epochs)):
            model.train()
            train_losses = []
            val_losses = []
            for batch in train_dataloader:
                x, y = batch
                x = x.to(config.device)
                y = y.to(config.device)
                optimizer.zero_grad()
                device = next(model.parameters()).device  
                output = model(x)
                if config.criterion == "regression":       loss = loss_fn(output.float(), y.float())
                elif config.criterion == "classification": loss = loss_fn(output.view(-1), y.float().view(-1))

                train_losses.append(loss.item())
                loss.backward()
                optimizer.step()
            model.eval()
            with torch.no_grad():
                for batch in val_dataloader:
                    x, y = batch
                    x = x.to(config.device)
                    y = y.to(config.device)
                    output = model(x)
                    loss = loss_fn(output.float(), y.float())                
                    val_losses.append(loss.item())
                epoch_train_losses.append(mean(train_losses))
                epoch_val_losses.append(mean(val_losses))
                if config.early_stopping:
                    early_stopping(mean(val_losses), model)
                    if early_stopping.early_stop:
                        print(f"Early stopping at epoch {epoch}")
                        early_stopping.load_best_model(model)
                        break
        
        logs["all_val_labels"].append(val_labels)
        logs["fold_train_losses"].append(epoch_train_losses)
        logs["fold_val_losses"].append(epoch_val_losses)
        train_actual, train_pred = get_predictions(model, train_dataloader, config)
        val_actual, val_pred = get_predictions(model, val_dataloader, config)
        logs["fold_actuals"].append(val_actual)
        logs["fold_preds"].append(val_pred)
        logs["fold_train_actuals"].append(train_actual)
        logs["fold_train_preds"].append(train_pred)
        if config.criterion == "regression": 
            correlations = {
                "train_spearman": compute_spearman_correlation(train_actual, train_pred) or 0,
                "val_spearman": compute_spearman_correlation(val_actual, val_pred) or 0,
                "train_pearson": compute_pearson_correlation(train_actual, train_pred) or 0,
                "val_pearson": compute_pearson_correlation(val_actual, val_pred) or 0
            }
            all_correlations.append(correlations)
            logs["all_correlation"].append(correlations)
        # plot_epoch_losses(epoch_train_losses, epoch_val_losses, filename, output_dir)
        # plot_fold_scatter(val_actual, val_pred, train_actual, train_pred, correlations, fold, output_dir)        

        losses["train_losses"].append(epoch_train_losses), losses["val_losses"].append(epoch_val_losses)
        fold_outputs["fold_train_actual"].append(train_actual), fold_outputs["fold_train_pred"].append(train_pred), fold_outputs["fold_val_actual"].append(val_actual), fold_outputs["fold_val_pred"].append(val_pred), 
        log_memory_usage()
        if config.save_model: torch.save(model.state_dict(), os.path.join(model_output_dir, f"model_fold_{fold}.pth"))
        # plot_correlations(logs["all_correlation"], output_dir)

    # Test set
    print("Evaluating test set")
    test_loader = DataLoader(test_dataset, batch_size=config.batch_size, shuffle=False)
    test_actual, test_pred = get_predictions(model, test_loader, config)

    if config.criterion == "regression":
        test_metrics = {
            "test_spearman": compute_spearman_correlation(test_actual, test_pred),
            "test_pearson": compute_pearson_correlation(test_actual, test_pred),
        }
        plot_all_folds_scatter(fold_outputs, all_correlations, config.output_dir, filename="all_folds_scatter")
        plot_test_scatter(test_actual, test_pred, config)
    elif config.criterion == "classification":
        test_metrics = {"test_confusion": confusion_matrix(test_actual, (torch.sigmoid(torch.tensor(test_pred)) > 0.5).int())}
        plot_confusion_matrix(fold_outputs, config, name=None)
        plot_test_confusion_matrix(test_actual, test_pred, config)


    logs["test_results"] = test_metrics

    # Save logs and fold_outputs as pkl
    plot_all_folds_losses(losses, output_dir, filename="all_folds")
    if config.criterion == "regression":     plot_all_folds_scatter(fold_outputs, all_correlations, output_dir, filename="all_folds_scatter")
    elif config.criterion == "classification": plot_confusion_matrix(fold_outputs, config, name=None)

    logs_file = os.path.join(output_dir, "logs.pkl")
    with open(logs_file, "wb") as f:
        pickle.dump(logs, f)
    print(f"Logs outputs saved to {logs_file}")
    fold_outputs_file = os.path.join(output_dir, "fold_outputs.pkl")
    with open(fold_outputs_file, "wb") as f:
        pickle.dump(fold_outputs, f)
    print(f"Fold outputs saved to {fold_outputs_file}")

    return model, logs, fold_outputs, losses

###############
### Options ###
############### 

# model input
one_hot_input = True
kmer_table_input = False
seq_features_input = False

# predicted value
label = "medium_GFP_fraction_cells_with_condensates"

# data split
balanced_split = True
group_clusters_split = False # split uses group kfold on mmseq similarity clusters
random_split = False   # randomly split parents and assign children into the same groups



##################
### RUN MODELS ###
##################


if balanced_split:
    ### Balanced Train Test Split -- Elena ###
    balanced_splits = pd.read_csv("/Users/clairehsieh/Library/CloudStorage/OneDrive-Personal/Documents/UCLA/rotations/kalli_kappel/model/kappel-lab-condensate-pred/simple_models/train_test_val_split_evaluation/9_19_25_base_split/codenSeq_data_with_base_split.csv")
    balanced_splits = balanced_splits.dropna(subset=[label])

    train = balanced_splits.loc[balanced_splits["base_split"] == "train"]
    test = balanced_splits.loc[balanced_splits["base_split"] == "test"]
    val = balanced_splits.loc[balanced_splits["base_split"] == "val"]


    if one_hot_input:
        ## Run Classification Model on One Hot Encoded Seq
        config = Config(output_dir = "/Users/clairehsieh/OneDrive/Documents/UCLA/rotations/kalli_kappel/results/logo/classification/balanced_split/",
                        batch_size=100, input_dim=1320, output_dim=1, criterion="classification")
        print(config)
        model_input = {"train":{'x':0, 'y':0}, 
                    "test":{'x':0, 'y':0}, 
                    "val":{'x':0, 'y':0}}
        for k, df in zip(model_input.keys(), [train, test, val]):
            model_input[k]['x'] = np.array([one_hot_encode(row['protein_seq']) for i,row in df.iterrows()])
            model_input[k]['x'] = torch.tensor(model_input[k]['x'].reshape(model_input[k]['x'].shape[0], -1), dtype=torch.float32).to(config.device)
            model_input[k]['y'] = (torch.tensor(np.array([row[label] for i, row in df.iterrows()]), dtype=torch.float32) >= 0.3).int().to(config.device)

        train_with_predefined_splits(model_input, config)

    if kmer_table_input:
        ### Run Classification Model on Kmer Tables ###
        config = Config(output_dir = "/Users/clairehsieh/OneDrive/Documents/UCLA/rotations/kalli_kappel/results/logo/classification/balanced_split/kmer/",
                        batch_size=100, input_dim=7696, output_dim=1, train_split="logo", criterion="classification")
        kmer_matrix, unique_kmers = kmer_table(list(balanced_splits["protein_seq"]), k=3)

        model_input = {"train":{'x':0, 'y':0}, 
                    "test":{'x':0, 'y':0}, 
                    "val":{'x':0, 'y':0}}
        for k, df, n in zip(model_input.keys(), [train, test, val], ["train", "test", "val"]):
            model_input[k]['y'] = (torch.tensor(np.array([row[label] for i, row in df.iterrows()]), dtype=torch.float32) >= 0.3).int().to(config.device)
            model_input[k]['x'] = torch.tensor(kmer_matrix[balanced_splits["base_split"] == n], dtype=torch.float32).to(config.device)

        # 6409 sequences, 7696 unique kmers with k = 3
        train_with_predefined_splits(model_input, config)

    if seq_features_input:
        ## Run Classification Model on Sequence Features ###
        config = Config(output_dir = "/Users/clairehsieh/OneDrive/Documents/UCLA/rotations/kalli_kappel/results/logo/classification/balanced_split/seq_features_aa_only/",
                        batch_size=100, input_dim=20, output_dim=1, criterion="classification")

        sequence_feature_cols = ['fraction_A',
            'fraction_C',
            'fraction_D',
            'fraction_E',
            'fraction_F',
            'fraction_G',
            'fraction_H',
            'fraction_I',
            'fraction_K',
            'fraction_L',
            'fraction_M',
            'fraction_N',
            'fraction_P',
            'fraction_Q',
            'fraction_R',
            'fraction_S',
            'fraction_T',
            'fraction_V',
            'fraction_W',
            'fraction_Y']
        model_input = {"train":{'x':0, 'y':0}, 
                    "test":{'x':0, 'y':0}, 
                    "val":{'x':0, 'y':0}}
        for k, df, n in zip(model_input.keys(), [train, test, val], ["train", "test", "val"]):
            model_input[k]['y'] = (torch.tensor(np.array([row[label] for i, row in df.iterrows()]), dtype=torch.float32) >= 0.3).int().to(config.device)
            model_input[k]['x'] = torch.tensor(df[sequence_feature_cols].values, dtype=torch.float32).to(config.device)

        train_with_predefined_splits(model_input, config)

if group_clusters_split:
    ## GROUPED BY CLUSTERS ###

    # Read in Train Test Split with Parent Child groups #
    # format: dict with protein_seq as id - match to group, label, x (features, one hot, or kmer table)

    clusters = pd.read_csv('/Users/clairehsieh/Library/CloudStorage/OneDrive-Personal/Documents/UCLA/rotations/kalli_kappel/data/parent_child_clusters.csv')
    cluster_ids = [int(k) for k, v in Counter(clusters["cluster_id"].values).items() if v > 2]
    clusters = clusters[clusters['cluster_id'].isin(cluster_ids)]

    if one_hot_input:
        model_input = {row["protein_seq"] : {"group":0, 
                                            "x":one_hot_encode(row["protein_seq"]), 
                                            "y":row[label]} for i,row in large_pl_data.dropna(subset=[label]).iterrows()}
    elif kmer_table_input:     
        kmer_matrix, unique_kmers = kmer_table(large_pl_data.dropna(subset=[label])["protein_seq"])
        tmp = {k:v for k,v in zip(large_pl_data.dropna(subset=[label])["protein_seq"], kmer_matrix)}
        model_input = {row["protein_seq"] : {"group":0, 
                                                "x" : tmp[row["protein_seq"]],
                                                "y " :row[label]} for i,row in large_pl_data.dropna(subset=[label]).iterrows()}
    elif seq_features_input: 
        sequence_feature_cols = ['fraction_A',
            'fraction_C',
            'fraction_D',
            'fraction_E',
            'fraction_F',
            'fraction_G',
            'fraction_H',
            'fraction_I',
            'fraction_K',
            'fraction_L',
            'fraction_M',
            'fraction_N',
            'fraction_P',
            'fraction_Q',
            'fraction_R',
            'fraction_S',
            'fraction_T',
            'fraction_V',
            'fraction_W',
            'fraction_Y',
            'ratio_R_K',
            'ratio_D_E',
            'ratio_S_G',
            'ratio_N_Q',
            'ratio_Y_F',
            'ratio_F_W',
            'ratio_Y_W',
            'ratio_R_Q',
            'ratio_K_Q',
            'fraction_group_ILMV',
            'fraction_group_RK',
            'fraction_group_DE',
            'fraction_group_GS',
            'fraction_group_YFW',
            'ratio_group_FYW_ILV',
            'ratio_group_FYW_R',
            'NCPR', # Net Charge Per Residue
            'FCR',
            'fraction_disorder_promoting',
            'mean_hydropathy',
            'frac_disorder',
            'avg_disorder',
            'sublibrary',
            'omega_STNQCH',
            'kappa_STNQCH_ILMV',
            'kappa_STNQCH_RK',
            'kappa_STNQCH_ED',
            'kappa_STNQCH_FWY',
            'kappa_STNQCH_A',
            'kappa_STNQCH_P',
            'kappa_STNQCH_G',
            'omega_ILMV',
            'kappa_ILMV_RK',
            'kappa_ILMV_ED',
            'kappa_ILMV_FWY',
            'kappa_ILMV_A',
            'kappa_ILMV_P',
            'kappa_ILMV_G',
            'omega_RK',
            'kappa_RK_ED',
            'kappa_RK_FWY',
            'kappa_RK_A',
            'kappa_RK_P',
            'kappa_RK_G',
            'omega_ED',
            'kappa_ED_FWY',
            'kappa_ED_A',
            'kappa_ED_P',
            'kappa_ED_G',
            'omega_FWY',
            'kappa_FWY_A',
            'kappa_FWY_P',
            'kappa_FWY_G',
            'omega_A',
            'kappa_A_P',
            'kappa_A_G',
            'omega_P',
            'kappa_P_G',
            'omega_G',
            'fraction_R_of_RK',
            'fraction_D_of_DE',
            'fraction_S_of_SG',
            'fraction_N_of_NQ',
            'fraction_Y_of_YF',
            'fraction_F_of_FW',
            'fraction_Y_of_YW',
            'fraction_R_of_RQ',
            'fraction_K_of_KQ',
            'fraction_group_ILV',
            'fraction_FYW_of_FYWILV',
            'fraction_FYW_of_FYWR']
        sequence_features = large_pl_data[["protein_seq", label] + sequence_feature_cols]
        sequence_features = sequence_features.set_index('protein_seq').loc[large_pl_data.dropna(subset=[label])["protein_seq"]].reset_index()
        model_input = {row["protein_seq"] : {"group":0, 
                                            "x" : sequence_features.loc[sequence_features["protein_seq"] == row['protein_seq']].drop(columns=["protein_seq", label]),
                                            "y" : row[label] } for i,row in large_pl_data.dropna(subset=[label]).iterrows() }

    for i, row in clusters.iterrows():
        try:        model_input[row["protein_seq"]]["group"] = row["cluster_id"]
        except:     pass

    if one_hot_input:
        ### Run Classification Model on One Hot Encoded Seq
        config = Config(output_dir = "/Users/clairehsieh/OneDrive/Documents/UCLA/rotations/kalli_kappel/results/logo/classification/group_parent_children_split/",
                        batch_size=10, input_dim=1320, output_dim=1, train_split="group_kfold", criterion="classification")
        print(config)
        x = torch.tensor([i["x"] for i in model_input.values()], dtype=torch.float32).to(config.device)
        x = x.reshape(x.shape[0], -1)
        y = (torch.tensor([i["y"] for i in model_input.values()], dtype=torch.float32) >= 0.3).int().to(config.device) # for classification
        group = [i["group"] for i in model_input.values()]

        model, logs, fold_outputs, losses = train(x, y, config, groups=group)

if random_split:
    ### Split parents randomly and assign children into the same groups ###
    # logo for folds
    dfs = [pd.read_csv(f"/Users/clairehsieh/Library/CloudStorage/OneDrive-Personal/Documents/UCLA/rotations/kalli_kappel/data/poolB_12ntBC_sublibrary{i}.csv") for i in range(1, 4)]

    parent_protein_seqs = (
        pd.concat([df[['base_protein_seq', 'protein_seq']] for df in dfs], ignore_index=True)
        .dropna()
    )
    # 13077 parent sequences

    base_seqs = (
        pd.concat([df[['base_protein_seq']] for df in dfs], ignore_index=True)
        .drop_duplicates(keep='first')
        .dropna()
    )
    # random split -- create dict in format {group : [list of parent seqs], ... }
    parent_seqs = list(base_seqs["base_protein_seq"])
    random.shuffle(parent_seqs)

    k = 10
    size = len(parent_seqs) // k
    groups = [parent_seqs[i*size:(i+1)*size] for i in range(k-1)]
    groups.append(parent_seqs[(k-1)*size:])  
    group_to_seq = {i:v for i,v in enumerate(groups)}

    # format into X, y, groups for logo split
    # groups should be group identity of each element in X (matched by index)
    subset_data = large_pl_data.merge(parent_protein_seqs, on="protein_seq")[["protein_seq", label, "base_protein_seq"]]
    subset_data = subset_data.dropna(subset=[subset_data.columns[1]])
    seq_data = np.array([one_hot_encode(i) for i in subset_data["protein_seq"]])
    labels = subset_data[label]
    # 6409

    # reverse group_to_seq dictionary --> {parent_seq1: group, ... } and map group onto child sequences
    seq_to_group = {s: g for g, lst in group_to_seq.items() for s in lst}
    subset_data["group"] = subset_data["base_protein_seq"].map(seq_to_group).tolist()
    subset_data = subset_data.dropna()
    groups_list = [seq_to_group[row["base_protein_seq"]] for ind, row in subset_data.iterrows()]

    if one_hot_input:
        config = Config(output_dir = "/Users/clairehsieh/OneDrive/Documents/UCLA/rotations/kalli_kappel/results/logo/classification/random_split/",
                        batch_size=10, input_dim=1320, output_dim=1, train_split="logo", criterion="classification")

        seq_data = torch.tensor(seq_data.reshape(seq_data.shape[0], -1), dtype=torch.float32).to(config.device)
        labels = (torch.tensor(labels, dtype=torch.float32) >= 0.5).int().to(config.device) # for classification

        model, logs, fold_outputs, losses = train(seq_data, labels, config, groups=groups_list)


    if kmer_table_input:
        seqs = subset_data["protein_seq"].values
        kmer_matrix, unique_kmers = kmer_table(seqs, k=3)
        # 6409 sequences, 7530 unique kmers with k = 3
        config = Config(output_dir = "/Users/clairehsieh/OneDrive/Documents/UCLA/rotations/kalli_kappel/results/logo/classification/random_split/kmer/",
                        batch_size=100, input_dim=7530, output_dim=1, train_split="logo", criterion="classification")
        kmer_matrix = torch.tensor(kmer_matrix, dtype=torch.float32).to(config.device)
        labels = (torch.tensor(labels, dtype=torch.float32) >= 0.3).int().to(config.device) # for classification
        model, logs, fold_outputs, losses = train(kmer_matrix, labels, config, groups=groups_list)

    
    if seq_features_input:
        ## Run Classification Model on Sequence Features ###
        config = Config(output_dir = "/Users/clairehsieh/OneDrive/Documents/UCLA/rotations/kalli_kappel/results/logo/classification/random_split/seq_features_aa_only/",
                        batch_size=100, input_dim=20, output_dim=1, criterion="classification")
        
        sequence_feature_cols = ['fraction_A',
            'fraction_C',
            'fraction_D',
            'fraction_E',
            'fraction_F',
            'fraction_G',
            'fraction_H',
            'fraction_I',
            'fraction_K',
            'fraction_L',
            'fraction_M',
            'fraction_N',
            'fraction_P',
            'fraction_Q',
            'fraction_R',
            'fraction_S',
            'fraction_T',
            'fraction_V',
            'fraction_W',
            'fraction_Y']

        x = torch.tensor(subset_data[sequence_feature_cols].values, dtype=torch.float32).to(config.device)
        labels = (torch.tensor(labels, dtype=torch.float32) >= 0.5).int().to(config.device) # for classification
        model, logs, fold_outputs, losses = train(x, labels, config, groups=groups_list)