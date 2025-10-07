#!/bin/bash
#$ -cwd
#$ -o /u/home/c/chsieh/job-output/joblog.$JOB_ID
#$ -j y
#$ -pe shared 2
#$ -l h_rt=8:00:00,h_data=4G
#$ -M $USER@mail
#$ -m bea


# load the job environment:
. /u/local/Modules/default/init/modules.sh
module use /u/project/CCN/apps/modulefiles

module load anaconda3
conda activate esm2

echo "Job started at: $(date)"
python3 /u/home/c/chsieh/rotation/model/get_esm.py
echo "Job finished at: $(date)"
