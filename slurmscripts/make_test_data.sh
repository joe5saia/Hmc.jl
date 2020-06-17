#!/bin/sh

julia --project=. /research/hmc/code/run_hmm.jl official 250 test
julia --project=. /research/hmc/code/run_hmm.jl official 251 test
julia --project=. /research/hmc/code/run_hmm.jl official 255 test


julia --project=. /research/hmc/code/run_hmm.jl alter 250 test
julia --project=. /research/hmc/code/run_hmm.jl alter 251 test
julia --project=. /research/hmc/code/run_hmm.jl alter 255 test