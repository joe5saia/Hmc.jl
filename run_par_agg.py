#!/bin/python3
import os

series = ["official"]
noise = [0.1, 1.0, 3.0]
siglen = [1, 12]
files = ["filtered_means", "filtered_state_probs", "filtered_variances", "filtered_trans_probs", "forecasts"]

for s in series:
    for n in noise:
        for l in siglen:
            path = '/moto/sscc/projects/biasedexpectations/data/output/signals_' + str(s) + '_noise_' + str(n) + '_len_' + str(l)
            for f in files:
                slurmcmd = 'srun -A sscc --mem-per-core 5G -c 1 -t 0-01:00:00 --output /moto/sscc/projects/biasedexpectaions/output/slurm/aggregate%J.txt'
                cmd = slurmcmd + " julia --project=/moto/sscc/projects/biasedexpectaions/ aggregate_par.jl " + path + " " + f
                os.system(cmd)

