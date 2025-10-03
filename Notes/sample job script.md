```
#!/bin/bash
#$ -cwd --> – Run the job in the current working directory (otherwise, default is `$HOME`).
# error = Merged with joblog 
#$ -o joblog.$JOB_ID --> – Redirect job output to a file named `joblog.<jobID>`.  
#$ -j y  --> merges stderr into stdout (so you only have one log file).
#$ -pe shared 2 --> – Request a parallel environment (`shared`) with **2 CPU cores**.
#$ -l h_rt=8:00:00,h_data=4G --> total memory = h_data * cores., 8 hours runtime and 4 GB RAM per core
#$ -M clhsieh84@g.ucla.edu
# Notify when
#$ -m bea

# load the job environment:
. /u/local/Modules/default/init/modules.sh
module use /u/project/CCN/apps/modulefiles

# Load the FSL module
module load anaconda3
conda activate esm2

# This is optional
# More info here: https://www.ccn.ucla.edu/wiki/index.php/Hoffman2:FSL 
export NO_FSL_JOBS=true

# Your script content goes here...
```