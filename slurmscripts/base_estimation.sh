#!/bin/sh

#SBATCH --account=sscc # The account name for the job.
#SBATCH --job-name=official_estimate # The job name.
#SBATCH --array=120-579
#SBATCH --time=6:00:00 # The time the job will take to run.
#SBATCH -n 1 # one cpu core
#SBATCH --mem=5G # The memory the job will use per cpu core.
#SBATCH --output=output/slurm/slurm_official_%a.txt
#SBATCH --chdir=/moto/sscc/projects/biasedexpectations
#SBATCH --mail-type=END
#SBATCH --mail-user=js4956@columbia.edu


module load julia

julia --project=. /moto/sscc/projects/biasedexpectations/code/run_hmm.jl official ${SLURM_ARRAY_TASK_ID}
