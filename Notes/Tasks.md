
# Sept 2025
## 9/19/25
- [x] simple MLP to predict `medium_GFP_fraction_cells_with_condensates`
- random train test split
- RESULT: ~0.78 correlation on validation set

## 9/23/25
- [x] MMSeqs2 to cluster sequences by similarity for train test split
	- `mmseqs module input_db output_db args [options]`
	- easy-cluster
		-  easy-cluster clusters entries from a FASTA/FASTQ file using the cascaded clustering algorithm.
		- `mmseqs easy-cluster examples/DB.fasta clusterRes tmp`
		- Here, examples/DB.fasta is the input file, clusterRes is the output, and tmp is for temporary
	- `mmseqs easy-cluster DB.fasta clusterRes tmp --min-seq-id 0.9 -c 0.9 --cov-mode 1`
		- since many sequences should be similar?
		- 
files
- [ ] GRFP draft
- [x] OR use parent seq to split (and make sure parent seqs aren't too similar)
- [ ] 
## 9/29/25
- [ ] research fluorescent imaging (FISH?)
- [ ] Go over onboarding doc
- [x] get keys
- [ ] Fill out presentation 
- [ ] read [paper ](https://www.sciencedirect.com/science/article/pii/S0092867418307311?via%3Dihub)
- [ ] ask Gerald and Ian for LORs for GRFP
- [ ] 
