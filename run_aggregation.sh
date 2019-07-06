#!/bin/sh
#SBATCH --account=sscc # The account name for the job.
#SBATCH --job-name=aggregate # The job name.
#SBATCH -c 1 # The number of cpus to use.
#SBATCH --time=1:00:00 # The time the job will take to run.
#SBATCH --mem=25G # The memory the job will use per cpu core.
#SBATCH --output=output/slurm/slurm_aggreagate.txt
#SBATCH --chdir=/moto/sscc/projects/biasedexpectations


module load julia

julia --project=. aggregate.jl official
julia --project=. aggregate.jl alter

