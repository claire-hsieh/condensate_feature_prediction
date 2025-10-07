```
#!/bin/bash
#$ -cwd
#$ -o /u/home/c/chsieh/job-output/joblog.$JOB_ID
#$ -j y
#$ -pe shared 2
#$ -l h_rt=8:00:00,h_data=4G
#$ -M $USER@mail
#$ -m bea


module load anaconda3
conda activate [ENV]

echo "Job started at: $(date)"

[SCRIPT CONTENT]


echo "Job finished at: $(date)"



```