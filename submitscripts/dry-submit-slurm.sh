#!/bin/bash -login
#SBATCH -p high                # partition, or queue, to assign to
#SBATCH -J 16s-dry             # name for job
#SBATCH -N 1                   # one "node", or computer
#SBATCH -n 1                   # one task for this node
#SBATCH -c 1                   # one core per task
#SBATCH -t 4:00:00             # 4 hours for run
#SBATCH --mem=20gb             # ask for no more than 10 GB of memory
#SBATCH --mail-user=nasiegel@ucdavis.edu        # email address to recieve notification
#SBATCH --mail-type=ALL                         # NOTIFICATIONS OF SLURM JOB STATUS

# initialize conda
. ~/miniconda3/etc/profile.d/conda.sh

# activate your desired conda environment
conda activate qiime2-snakemake

# fail on weird errors
set -o nounset
set -o errexit
set -x

# run the snakemake!
# Select which snakefile you want to submit
cd /home/nasiegel/2022-tagseq-qiime2-snakemake-2
snakemake -p -j 8 --use-conda --until multiqc

# print out various information about the job
env | grep SLURM            # Print out values of the current jobs SLURM environment variables

scontrol show job ${SLURM_JOB_ID}     # Print out final statistics about resource uses before job exits

sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch
