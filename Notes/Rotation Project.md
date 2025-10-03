## hoffmann
ssh chsieh@hoffman2.idre.ucla.edu                                   
FarewellHope93%
scp chsieh@hoffman2.idre.ucla.edu:/u/project/kappel/chsieh/rotation/


- https://www.hoffman2.idre.ucla.edu/Policies/Job-scheduling-policy.html
- [interactive session](https://www.ccn.ucla.edu/wiki/index.php/Hoffman2:Interactive_Sessions): `qrsh`
- run jobs
	- `qsub -cwd -V -N J1 -l h_data=64M,h_rt=00:05:00 -M $USER -m bea /u/project/CCN/apps/examples/qsub/gather.sh`
	- interactive
		- qlogin -l h_data=4G,h_rt=2:00:00
- [[sample job script]]
- to submit: `qsub job.sh`
- check status: `qstat -u chsieh`
- delete job: `qdel <jobID>`

# Data
- `sublibrary` = batch
- `natural_protein_seq_set, composition_seq_set, patterning_seq_set` -> where data came from
	- evenly split this
- First, a set of 79 “base sequences” was selected for extensive mutagenesis. These sequences include 43 of the fragments of natural protein sequences and one designed variant from the small sequence library, for a total of 44 “Class 1 base sequences.” Another 35 fragments of natural protein sequences (“Class 2 base sequences”) were further selected from PhasePro, excluding protein sequences that were annotated as partner dependent. To this end, all protein regions annotated to drive phase separation were examined, and all possible 66 amino acid fragments with amino acid composition and dipeptide composition similarity (Pearson correlation coefficient (r2)) to the Class 1 base sequences and to each other of less than 0.6 (“Class 3 sequences”) were identified. 35 of these sequences, including at most one sequence fragment per protein and prioritizing sequence fragments with the highest amino acid composition similarity to the full protein region as well as regions that were predicted to be more disordered.
- The natural protein sequence fragment set includes: (1) all base sequences, (2) all remaining Class 3 sequences (283 sequences), (3) Class 4A sequences (519 sequences): fragments from disordered (as annotated by MobiDB36) sequences from LLPSDB37 that were annotated either as phase separating or not phase separating with maximum amino acid composition and dipeptide correlation of 0.8 to each other and to all base sequences and all Class 3 sequences; and (4) Class 4B sequences (798 sequences): disordered regions (as annotated by MobiDB) from Disprot38 (release 2022_03) with maximum amino acid composition and dipeptide correlation of 0.6 to each other and to all base sequences, Class 3 sequences, and Class 4A sequences.
- large dataset columns
	- sequence
		- `'protein_seq', 'DNA_seq_1', 'DNA_seq_2', 'DNA_seq_3', 'DNA_seq_4', 'DNA_seq_5'`
		- `'barcode_8_1', 'barcode_8_2', 'barcode_8_3', 'barcode_8_4', 'barcode_8_5', 'sublibrary_tag_1'`
		- `'sublibrary_tag_2', 'sublibrary_tag_3', 'sublibrary_tag_4', 'sublibrary_tag_5', 'natural_protein_seq_set', 'composition_seq_set', 'patterning_seq_set'`
	- GFP : 
		- `'low_GFP_cell_counts', 'low_GFP_num_condensates', 'low_GFP_fraction_gfp_in_condensates',`
		- `'low_GFP_total_condensate_area', 'low_GFP_mean_condensate_area', 'low_GFP_std_condensate_area',`
		- `'low_GFP_mean_condensate_eccentricity', 'low_GFP_std_condensate_eccentricity',`
		- `'low_GFP_fraction_cells_with_condensates', 'low_GFP_glcm_dissimilarity', 'low_GFP_correlation_GFP_dapi',`
		- `'low_GFP_fraction_nucleolar', 'low_GFP_fraction_chromatin',`
		- `'medium_GFP_cell_counts', 'medium_GFP_num_condensates', 'medium_GFP_fraction_gfp_in_condensates', 'medium_GFP_total_condensate_area', 'medium_GFP_mean_condensate_area', 'medium_GFP_std_condensate_area', 'medium_GFP_mean_condensate_eccentricity', 'medium_GFP_std_condensate_eccentricity', 'medium_GFP_fraction_cells_with_condensates', 'medium_GFP_glcm_dissimilarity', 'medium_GFP_correlation_GFP_dapi', 'medium_GFP_fraction_nucleolar', 'medium_GFP_fraction_chromatin', 'high_GFP_cell_counts', 'high_GFP_num_condensates', 'high_GFP_fraction_gfp_in_condensates', 'high_GFP_total_condensate_area', 'high_GFP_mean_condensate_area', 'high_GFP_std_condensate_area', 'high_GFP_mean_condensate_eccentricity', 'high_GFP_std_condensate_eccentricity', 'high_GFP_fraction_cells_with_condensates', 'high_GFP_glcm_dissimilarity', 'high_GFP_correlation_GFP_dapi', 'high_GFP_fraction_nucleolar', 'high_GFP_fraction_chromatin',`
	- SNAP:
		- `'low_SNAP_cell_counts', 'low_SNAP_num_condensates', 'low_SNAP_fraction_gfp_in_condensates', 'low_SNAP_total_condensate_area', 'low_SNAP_mean_condensate_area', 'low_SNAP_std_condensate_area', 'low_SNAP_mean_condensate_eccentricity', 'low_SNAP_std_condensate_eccentricity', 'low_SNAP_fraction_cells_with_condensates', 'low_SNAP_glcm_dissimilarity', 'low_SNAP_correlation_GFP_dapi', 'low_SNAP_fraction_nucleolar', 'low_SNAP_fraction_chromatin', 'medium_SNAP_cell_counts', 'medium_SNAP_num_condensates', 'medium_SNAP_fraction_gfp_in_condensates', 'medium_SNAP_total_condensate_area', 'medium_SNAP_mean_condensate_area', 'medium_SNAP_std_condensate_area', 'medium_SNAP_mean_condensate_eccentricity', 'medium_SNAP_std_condensate_eccentricity', 'medium_SNAP_fraction_cells_with_condensates', 'medium_SNAP_glcm_dissimilarity', 'medium_SNAP_correlation_GFP_dapi', 'medium_SNAP_fraction_nucleolar', 'medium_SNAP_fraction_chromatin', 'high_SNAP_cell_counts', 'high_SNAP_num_condensates', 'high_SNAP_fraction_gfp_in_condensates', 'high_SNAP_total_condensate_area', 'high_SNAP_mean_condensate_area', 'high_SNAP_std_condensate_area', 'high_SNAP_mean_condensate_eccentricity', 'high_SNAP_std_condensate_eccentricity', 'high_SNAP_fraction_cells_with_condensates', 'high_SNAP_glcm_dissimilarity', 'high_SNAP_correlation_GFP_dapi', 'high_SNAP_fraction_nucleolar', 'high_SNAP_fraction_chromatin',`
		- SNAP:  protein that, when fused to another protein of interest, allows for specific and irreversible covalent labeling with a variety of synthetic probes, such as fluorescent dyes or biotin
	- `'base_sequence_gene_name',`
	-  amino acids
		- `'fraction_A', 'fraction_C', 'fraction_D', 'fraction_E', 'fraction_F', 'fraction_G', 'fraction_H', 'fraction_I', 'fraction_K', 'fraction_L', 'fraction_M', 'fraction_N', 'fraction_P', 'fraction_Q', 'fraction_R', 'fraction_S', 'fraction_T', 'fraction_V', 'fraction_W', 'fraction_Y', 'ratio_R_K', 'ratio_D_E', 'ratio_S_G', 'ratio_N_Q', 'ratio_Y_F', 'ratio_F_W', 'ratio_Y_W', 'ratio_R_Q', 'ratio_K_Q', 'fraction_group_ILMV', 'fraction_group_RK', 'fraction_group_DE', 'fraction_group_GS', 'fraction_group_YFW', 'ratio_group_FYW_ILV', 'ratio_group_FYW_R', 'NCPR', 'FCR', 'fraction_disorder_promoting', 'mean_hydropathy', 'frac_disorder', 'avg_disorder', 'sublibrary', 'omega_STNQCH', 'kappa_STNQCH_ILMV', 'kappa_STNQCH_RK', 'kappa_STNQCH_ED', 'kappa_STNQCH_FWY', 'kappa_STNQCH_A', 'kappa_STNQCH_P', 'kappa_STNQCH_G', 'omega_ILMV', 'kappa_ILMV_RK', 'kappa_ILMV_ED', 'kappa_ILMV_FWY', 'kappa_ILMV_A', 'kappa_ILMV_P', 'kappa_ILMV_G', 'omega_RK', 'kappa_RK_ED', 'kappa_RK_FWY', 'kappa_RK_A', 'kappa_RK_P', 'kappa_RK_G', 'omega_ED', 'kappa_ED_FWY', 'kappa_ED_A', 'kappa_ED_P', 'kappa_ED_G', 'omega_FWY', 'kappa_FWY_A', 'kappa_FWY_P', 'kappa_FWY_G', 'omega_A', 'kappa_A_P', 'kappa_A_G', 'omega_P', 'kappa_P_G', 'omega_G', 'fraction_R_of_RK', 'fraction_D_of_DE', 'fraction_S_of_SG', 'fraction_N_of_NQ', 'fraction_Y_of_YF', 'fraction_F_of_FW', 'fraction_Y_of_YW', 'fraction_R_of_RQ', 'fraction_K_of_KQ', 'fraction_group_ILV', 'fraction_FYW_of_FYWILV', 'fraction_FYW_of_FYWR']`
	- terms
		- `glcm` : gray level cooccurrence matrix
		- `dapi`: fluorescent nuclear stain 4′,6-diamidino-2-phenylindole., nuclear marker that provides ground-truth structure for separating nuclei in the images.
		- `low, medium, high` : concentration bins?
		- `eccentricity`: how elongated or irregular a condensate is relative to a perfect circle/sphere (0 = perfect circle)
		- `NCPR`:  (Net Charge Per Residue):
	- Distribution of `medium_GFP_fraction_cells_with_condensates`
		- ![[Pasted image 20251002083206.png|550]]
# Results
## 9/19/25
- naive results with random train test split for predicting `medium_GFP_fraction_cells_with_condensates`
	- ![[Pasted image 20250920101904.png|1000]]
	- ![[Pasted image 20250920101848.png|1000]]
## 9/23/25
### mmseqs all by all seq similarity
```bash
mmseqs createdb ../data/parent_seq_1.fasta ../data/mmseqs/parent_seq_1_DB
mmseqs search ../data/mmseqs/parent_seq_1_DB ../data/mmseqs/parent_seq_1_DB ../data/mmseqs/parent_seq__1_resultDB tmp --min-seq-id 0.3 -c 0.5
mmseqs convertalis ../data/mmseqs/parent_seq_1_DB ../data/mmseqs/parent_seq_1_DB ../data/mmseqs/parent_seq__1_resultDB  ../data/mmseqs/parent_seq__1_result.tsv
```


### clustering parent seqs 

| Hierarchical                              | DBScan                                    |
| ----------------------------------------- | ----------------------------------------- |
| ![[Pasted image 20250923105254.png\|550]] | ![[Pasted image 20250924081133.png\|550]] |

- using seq similarity scores from mmseqs2 as distance matrix, 
- **hierarchical**
	- 320 parent seqs
	- if i set number of clusters to be low (~10 or 20), most of the sequences will end up in one cluster
		- `Counter({10: 264, 1: 12, 2: 12, 5: 8, 9: 4, 4: 4, 3: 4, 7: 4, 8: 4, 6: 4})`
	- threshold = 0.9 --> 87 clusters, most with 4 members		
-  **DBscan**
	- `Counter({-1: 288, 0: 12, 1: 12, 2: 8})`
- parent seq clustering --> imbalanced ?? 
### clustering all sequences

| Hierarchical                              | DBScan                                    |
| ----------------------------------------- | ----------------------------------------- |
| ![[Pasted image 20250923111541.png\|550]] | ![[Pasted image 20250924150923.png\|550]] |

- **Hierarchical**
	- min-seq-id 0.3 -c 0.5
	- `Counter({10:12947, 1:367, 2:336, 8:149, 3:143, 6:140, 4:139, 9:139, 5:133, 7:129})`
- **DBScan**s
	- `Counter({-1: 5283, 50: 367, 2: 337, 11: 252, 36: 149, 18: 143, 16: 140, 43: 140, 1: 139,`
- random split with parent sequences as groups --> 10 groups each with 32 parent sequences
### logo train test split
- ![[Pasted image 20250923211129.png|5000]]
- ![[Pasted image 20250923211148.png|4000]]


## 10/2/25
- Classfication : random parent split, one hot encoding of sequence
	- ![[Pasted image 20251002090127.png|4000]]

### group k fold with parent seq clusters (from mia)
- one hot encoding
![[Pasted image 20251002134607.png|2000]]

### Balanced dataset (elena)

![[Pasted image 20251002211511.png|1000]]
![[Pasted image 20251002211520.png|2000]]![[Pasted image 20251002212301.png|550]]