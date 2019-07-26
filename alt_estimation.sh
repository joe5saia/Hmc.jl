#!/bin/sh

#SBATCH --account=sscc # The account name for the job.
#SBATCH --job-name=alt_estimate # The job name.
#SBATCH --array=121-578
#SBATCH --time=4:00:00 # The time the job will take to run.
#SBATCH -n 1 # one cpu core
#SBATCH --mem=5G # The memory the job will use per cpu core.
#SBATCH --output=output/slurm/slurm_alt_%a.txt
#SBATCH --chdir=/moto/sscc/projects/biasedexpectations


module load julia

julia --project=. /moto/sscc/projects/biasedexpectations/code/run_hmm.jl alter ${SLURM_ARRAY_TASK_ID}
