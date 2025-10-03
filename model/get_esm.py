import pathlib
import torch

from esm import FastaBatchedDataset, pretrained

# to install esm2
# pip install fair-esm  # latest release, OR:
# pip install git+https://github.com/facebookresearch/esm.git  # bleeding edge, current repo main branch

# Checkpoint name	Num layers	Num parameters
# esm2_t48_15B_UR50D	48	15B
# esm2_t36_3B_UR50D	    36	3B
# esm2_t33_650M_UR50D	33	650M
# esm2_t30_150M_UR50D	30	150M
# esm2_t12_35M_UR50D	12	35M
# esm2_t6_8M_UR50D	    6	8M

# from kaggle: https://www.kaggle.com/code/viktorfairuschin/extracting-esm-2-embeddings-from-fasta-files
def extract_embeddings(model_name, fasta_file, output_dir, tokens_per_batch=4096, seq_length=1022,repr_layers=[33]):
    
    model, alphabet = pretrained.load_model_and_alphabet(model_name)
    model.eval()

    if torch.cuda.is_available():
        model = model.cuda()
        
    dataset = FastaBatchedDataset.from_file(fasta_file)
    batches = dataset.get_batch_indices(tokens_per_batch, extra_toks_per_seq=1)

    data_loader = torch.utils.data.DataLoader(
        dataset, 
        collate_fn=alphabet.get_batch_converter(seq_length), 
        batch_sampler=batches
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    
    with torch.no_grad():
        for batch_idx, (labels, strs, toks) in enumerate(data_loader):

            print(f'Processing batch {batch_idx + 1} of {len(batches)}')

            if torch.cuda.is_available():
                toks = toks.to(device="cuda", non_blocking=True)

            out = model(toks, repr_layers=repr_layers, return_contacts=False)

            logits = out["logits"].to(device="cpu")
            representations = {layer: t.to(device="cpu") for layer, t in out["representations"].items()}
            
            for i, label in enumerate(labels):
                entry_id = label.split()[0]
                
                filename = output_dir / f"{entry_id}.pt"
                truncate_len = min(seq_length, len(strs[i]))

                result = {"entry_id": entry_id}
                result["mean_representations"] = {
                        layer: t[i, 1 : truncate_len + 1].mean(0).clone()
                        for layer, t in representations.items()
                    }

                torch.save(result, filename)

model_name = 'esm2_t6_8M_UR50D'
fasta_file = pathlib.Path('/u/project/kappel/chsieh/rotation/data/parent_seq.fasta')
output_dir = pathlib.Path('/u/project/kappel/chsieh/rotation/data/esm2/parent/')

extract_embeddings(model_name, fasta_file, output_dir)

