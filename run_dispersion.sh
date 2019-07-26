#!/bin/sh
#SBATCH --account=sscc # The account name for the job.
#SBATCH --job-name=dispersion # The job name.
#SBATCH -c 1 # The number of cpus to use.
#SBATCH --time=2:00:00 # The time the job will take to run.
#SBATCH --mem=10G # The memory the job will use per cpu core.
#SBATCH --chdir=/moto/sscc/projects/biasedexpectations
#SBATCH --output=output/slurm/slurm_dispersion.txt


module load julia

julia --project=. code/dispersion.jl

