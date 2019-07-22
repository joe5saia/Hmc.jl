#!/bin/sh

julia --project=. /research/hmc/code/run_hmm.jl official 250
julia --project=. /research/hmc/code/run_hmm.jl official 251
julia --project=. /research/hmc/code/run_hmm.jl official 255


julia --project=. /research/hmc/code/run_hmm.jl alter 250
julia --project=. /research/hmc/code/run_hmm.jl alter 251
julia --project=. /research/hmc/code/run_hmm.jl alter 255