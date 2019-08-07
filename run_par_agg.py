#!/bin/python3
import os

series = ["official"]
noise = [0.1, 1.0, 3.0]
siglen = [1, 12]
files = ["filtered_means", "filtered_state_probs", "filtered_variances", "filtered_trans_probs", "forecasts"]

for s in series:
    for n in noise:
        for l in siglen:
            #path = '/moto/sscc/projects/biasedexpectations/data/output/signals_' + str(s) + '_noise_' + str(n) + '_len_' + str(l)
            path = '/moto/sscc/projects/biasedexpectations/data/output/signals_' + str(s) + '_noise_' + str(n) + '_allsignal'
            for f in files:
                slurmcmd = 'srun -A sscc --mem-per-cpu 5G -c 1 -t 0-01:00:00 --output /moto/sscc/projects/biasedexpectations/output/slurm/aggregate%J.txt'
                cmd = slurmcmd + " julia --project=/moto/sscc/projects/biasedexpectations/ /moto/sscc/projects/biasedexpectations/code/aggregate_par.jl " + path + " " + f + " &"
                os.system(cmd)

