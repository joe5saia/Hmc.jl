#!/bin/sh

#SBATCH --account=sscc # The account name for the job.
#SBATCH --job-name=insample # The job name.
#SBATCH --time=4:00:00 # The time the job will take to run.
#SBATCH -n 1 # one cpu core
#SBATCH --mem=5G # The memory the job will use per cpu core.
#SBATCH --output=output/slurm/slurm_insample.txt
#SBATCH --chdir=/moto/sscc/projects/biasedexpectations
#SBATCH --mail-type=END
#SBATCH --mail-user=js4956@columbia.edu


module load julia

julia --project=. /moto/sscc/projects/biasedexpectations/code/run_insamplefcasts.jl 
