#!/bin/bash -login
#SBATCH -p high                # partition, or queue, to assign to
#SBATCH -J mld                 # name for job
#SBATCH -N 1                   # one "node", or computer
#SBATCH -n 1                   # one task for this node
#SBATCH -c 1                   # one core per task
#SBATCH -t 24:00:00             # ask for no more than 30 minutes
#SBATCH --mem=20Gb             # 20Gb should be enough

# initialize conda
. ~/miniconda3/etc/profile.d/conda.sh

# activate your desired conda environment
conda activate snake-tagseq

# go to the directory you ran 'sbatch' in, OR just hardcode it...

# fail on weird errors
set -o nounset
set -o errexit
set -x

# run the snakemake!
# Select which snakefile you want to submit
cd /home/nasiegel/MLD/tagseq-qiime2-snakemake-1
snakemake --use-conda

# print out various information about the job
env | grep SLURM            # Print out values of the current jobs SLURM environment variables

scontrol show job ${SLURM_JOB_ID}     # Print out final statistics about resource uses before job exits

sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch
