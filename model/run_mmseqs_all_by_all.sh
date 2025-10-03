#!/bin/bash
# Usage: ./script.sh tmp_dir [min_seq_id] [coverage] [extra mmseqs args]

tmp_dir=$1
input_file=$2
min_seq_id=${3:-0.3}   # default 0.3
coverage=${4:-0.8}     # default 0.8
shift 4 || true        # drop first three args if present
extra_args="$@"        # capture remaining args

mmseqs createdb $input_file "$tmp_dir/tmp_DB"

mmseqs search "$tmp_dir/tmp_DB" "$tmp_dir/tmp_DB" "$tmp_dir/tmp_result_DB" tmp \
    --min-seq-id "$min_seq_id" -c "$coverage" $extra_args

mmseqs convertalis "$tmp_dir/tmp_DB" "$tmp_dir/tmp_DB" "$tmp_dir/tmp_result_DB" "$tmp_dir/tmp_result.tsv"


# result file headers: 
# Query identifier 	Target identifier 	Sequence identity 	Alignment length 	Number of mismatches 	Number of gap openings 	Query start position 	Query end position 	Target start position 	Target end position 	E-value 	Bit score
