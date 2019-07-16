#!/bin/sh

#SBATCH --account=sscc # The account name for the job.
#SBATCH --job-name=alt_estimate # The job name.
#SBATCH --array=121-546 
#SBATCH --time=1-0:00:00 # The time the job will take to run.
#SBATCH -n 1 # one cpu core
#SBATCH --mem=20G # The memory the job will use per cpu core.
#SBATCH --output=output/slurm/slurm_official_%a.txt
#SBATCH --chdir=/moto/sscc/projects/biasedexpectations


module load julia

julia --project=. /moto/sscc/projects/biasedexpectations/code/run_hmm.jl alter %a
