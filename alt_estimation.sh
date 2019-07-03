#!/bin/sh
#
#SBATCH --account=sscc # The account name for the job.
#SBATCH --job-name=runestiamtion_alt # The job name.
#SBATCH -c 1 # The number of cpu cores to use.
#SBATCH --time=10:00:00 # The time the job will take to run.
#SBATCH --mem-per-cpu=5gb # The memory the job will use per cpu core.
#SBATCH --output=output/slurm/slurm_alt.txt
#SBATCH --chdir=/moto/sscc/projects/biasedexpectations


module load julia

julia --project=. run_hmm.jl alter
 

#End of script
