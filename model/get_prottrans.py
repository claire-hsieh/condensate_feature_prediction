from transformers import T5Tokenizer, T5EncoderModel
import torch
import re
import os
import time

"""
Available models

ProtT5 variants
    Rostlab/prot_t5_xl_uniref50 → full XL model
    Rostlab/prot_t5_xl_half_uniref50-enc → smaller XL “half” variant (smaller memory footprint)
    Rostlab/prot_t5_base_mt_uniref50 → base model, much smaller than XL
    Rostlab/prot_t5_xxl_uniref50 → XXL model
    Rostlab/prot_t5_xl_bfd → XL trained on BFD
    Rostlab/prot_t5_xxl_bfd → XXL trained on BFD
ProtBert variants
    Rostlab/prot_bert
    Rostlab/prot_bert_bfd
    Rostlab/prot_bert_bfd_ss3
    Rostlab/prot_bert_bfd_membrane
    Rostlab/prot_bert_bfd_localization
Other transformer models
    Rostlab/prot_albert
    Rostlab/prot_xlnet
    Rostlab/prot_electra_generator_bfd
    Rostlab/prot_electra_discriminator_bfd
"""


def extract_prottrans_embeddings(model_name, fasta_file, output_dir, device='cuda'):
    """
    fasta_file: path to FASTA
    output_dir: Path object
    """
    # Load model and tokenizer
    tokenizer = T5Tokenizer.from_pretrained(model_name, do_lower_case=False)
    model = T5EncoderModel.from_pretrained(model_name, use_safetensors=True)

    model.eval()
    if device == 'cuda' and torch.cuda.is_available():
        model = model.to('cuda')

    # Read FASTA sequences
    sequences = {}
    with open(fasta_file) as f:
        current_id = None
        current_seq = []
        for line in f:
            if line.startswith(">"):
                if current_id:
                    sequences[current_id] = "".join(current_seq)
                current_id = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_id:
            sequences[current_id] = "".join(current_seq)

    # Process sequences one by one (can be batched)
    with torch.no_grad():
        for seq_id, seq in sequences.items():
            # ProtTrans expects spaces between amino acids
            seq_spaced = " ".join(list(seq))
            inputs = tokenizer(seq_spaced, return_tensors="pt")
            if device == 'cuda' and torch.cuda.is_available():
                inputs = {k:v.to('cuda') for k,v in inputs.items()}

            outputs = model(**inputs)
            # Last hidden state: (batch, seq_len, hidden_dim)
            embedding = outputs.last_hidden_state.squeeze(0).cpu()
            
            torch.save({"entry_id": seq_id, "embedding": embedding}, f"{output_dir}/{seq_id}.pt")
            print(f"Saved embedding for {seq_id}")


time1 = time.time()
fasta_file = "/Users/clairehsieh/Library/CloudStorage/OneDrive-Personal/Documents/UCLA/rotations/kalli_kappel/data/parent_seqs.fasta"
output_dir = "/Users/clairehsieh/Library/CloudStorage/OneDrive-Personal/Documents/UCLA/rotations/kalli_kappel/data/prottrans/parent/"

# fasta_file = "/u/home/c/chsieh/rotation/data/test.fasta"
# output_dir = "/u/home/c/chsieh/rotation/data/prottrans/"

if not os.path.exists(output_dir):
    os.makedirs(output_dir, exist_ok=True)

model_name = "Rostlab/prot_t5_base_mt_uniref50"  # or prot_t5_base_mt_uniref50 for smaller model

extract_prottrans_embeddings(model_name, fasta_file, output_dir)

time2 = time.time()
print(f"Time taken: {time2 - time1}")




# device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

# # Load the tokenizer
# tokenizer = T5Tokenizer.from_pretrained('Rostlab/prot_t5_base_mt_uniref50', do_lower_case=False)

# # Load the model
# model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_base_mt_uniref50").to(device)

# # only GPUs support half-precision currently; if you want to run on CPU use full-precision (not recommended, much slower)
# if device == torch.device("cpu"):
#     model.to(torch.float32)

# # prepare your protein sequences as a list
# sequence_examples = ["PRTEINO", "SEQWENCE"]

# # replace all rare/ambiguous amino acids by X and introduce white-space between all amino acids
# sequence_examples = [" ".join(list(re.sub(r"[UZOB]", "X", sequence))) for sequence in sequence_examples]

# # tokenize sequences and pad up to the longest sequence in the batch
# ids = tokenizer(sequence_examples, add_special_tokens=True, padding="longest")

# input_ids = torch.tensor(ids['input_ids']).to(device)
# attention_mask = torch.tensor(ids['attention_mask']).to(device)

# # generate embeddings
# with torch.no_grad():
#     embedding_repr = model(input_ids=input_ids, attention_mask=attention_mask)

# # extract residue embeddings for the first ([0,:]) sequence in the batch and remove padded & special tokens ([0,:7]) 
# emb_0 = embedding_repr.last_hidden_state[0,:7] # shape (7 x 1024)
# # same for the second ([1,:]) sequence but taking into account different sequence lengths ([1,:8])
# emb_1 = embedding_repr.last_hidden_state[1,:8] # shape (8 x 1024)

# # if you want to derive a single representation (per-protein embedding) for the whole protein
# emb_0_per_protein = emb_0.mean(dim=0) # shape (1024)