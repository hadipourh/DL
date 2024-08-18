#!/bin/bash

#SBATCH --nodelist=xeon192g0                          # Host's name (request a specific list of hosts)
#SBATCH --job-name=warp                          # Name of job
#SBATCH --mail-type=ALL                              # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=hossein.hadipour@iaik.tugraz.at     # Where to send mail	
# SBATCH --ntasks=1                                      # Run on a single CPU (number of tasks to run)
# SBATCH --mem=120gb                                      # Job memory request (minimum amount of real memory)
# SBATCH --time=72:59:59                                 # Time limit hrs:min:sec
#SBATCH --output=/workstore/hhadipour/boomerang/warp_%A_%a.log                     # Standard output log
#SBATCH --error=/workstore/hhadipour/boomerang/warp_%A_%a.err                     # Standard error log
#SBATCH --cpus-per-task=1                              # Number of cpus required per task
pwd; hostname; date
# Install python-sat
# make
./difflin ${SLURM_ARRAY_TASK_ID}
date

## sbatch --array=0-63 -N1 jobfile.sh

