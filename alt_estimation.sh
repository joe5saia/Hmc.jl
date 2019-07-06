#!/bin/sh
#
#SBATCH --account=sscc # The account name for the job.
#SBATCH --job-name=alt_estimate # The job name.
#SBATCH -N 1 # The number of nodes to use.
#SBATCH --time=05-0:00:00 # The time the job will take to run.
#SBATCH --mem=100G # The memory the job will use per cpu core.
#SBATCH --output=output/slurm/slurm_alt.txt
#SBATCH --chdir=/moto/sscc/projects/biasedexpectations


module load julia

julia --project=. -p 24 run_hmm.jl alter
julia --project=. aggregate.jl alter
