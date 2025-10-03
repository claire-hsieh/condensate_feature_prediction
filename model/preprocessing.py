import json
import pandas as pd

### Parent Child Seqs Processing ###
# fasta_id_seq.json -- dictionary matching parent and child fasta id to seq; since protien ids are duplicated, i just used numbers, but fasta_id_seq[parent][1] is not the same as fasta_id_seq[child][1]
# parent_seq.fasta -- fasta of all parent seqs (from poolB... files)
# child_seq.fasta -- fasta of all child seqs (from large_pool_data... )

large_pl_data = pd.read_csv("/Users/clairehsieh/OneDrive/Documents/UCLA/rotations/kalli_kappel/data/large_pool_data_with_seq_info_202507.csv")
dfs = [pd.read_csv(f"/Users/clairehsieh/Library/CloudStorage/OneDrive-Personal/Documents/UCLA/rotations/kalli_kappel/data/poolB_12ntBC_sublibrary{i}.csv") for i in range(1, 4)]
base_seqs = (
    pd.concat([df[['protein_id', 'base_protein_seq']] for df in dfs], ignore_index=True)
      .drop_duplicates(keep='first')
      .dropna()
)

base_seqs = base_seqs.reset_index(drop=True)
base_seqs = base_seqs[~base_seqs['protein_id'].duplicated(keep='first')]

fasta_id_seq = {"parent":{}, 
                "child": {}}

with open("/Users/clairehsieh/Library/CloudStorage/OneDrive-Personal/Documents/UCLA/rotations/kalli_kappel/data/child_seqs.fasta", "w") as w:
    for i, row in large_pl_data.iterrows(): 
        w.write(f">{i}\n")
        w.write(f"{row["protein_seq"]}\n")
        fasta_id_seq['child'][i] = row['protein_seq']

with open("/Users/clairehsieh/Library/CloudStorage/OneDrive-Personal/Documents/UCLA/rotations/kalli_kappel/data/fasta_id_seq.json", "w") as f:
    json.dump(fasta_id_seq, f)

