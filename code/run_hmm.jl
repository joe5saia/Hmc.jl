###############################################################################
#= 
Program: run_hmm.jl 
Author: Joe Saia
Date: July 2019
Info: Runs code to estimate hidden markov models with Gibbs sampling as 
defined in the module hmc.jl for a range of horizions and saves results
Input: Takes a Stata file with the inflation series as an Input
Output: Writes several csvs
=#
###############################################################################
using Distributed

## CHANGE THIS TO DIRECTORY LOCATION
@everywhere root_dir = "/research/hmc/"
cd(root_dir)

## Supply the data series as a command line argument 
if length(ARGS) == 0
    @error "Must supply either \"official\" or \"alter\" as command line argument e.g. \n\tjulia --project code/run_hmm.jl official"
    exit(1)
elseif ARGS[1] == "official"
    series = :offic_inf
    filesuffix = "official"
    println("Using official inflation")
elseif ARGS[1] == "alter"
    series = :alter_inf
    filesuffix = "alter"
    println("Using alternative inflation")
else 
    @error "Must supply either \"official\" or \"alter\" as command line argument e.g. \n\tjulia --project code/run_hmm.jl official"
    exit(2)
end


@everywhere using Pkg
@everywhere Pkg.activate(root_dir)
@everywhere push!(LOAD_PATH, "$(root_dir)src")
@everywhere using Hmc

using StatFiles
using DataFrames
using Dates


## Make sure that paths exist and if not create them
!ispath("data/output/$(filesuffix)/") && mkpath("data/output/$(filesuffix)/")
!ispath("data/output/signals/$(filesuffix)/low") && mkpath("data/output/signals/$(filesuffix)/low")
!ispath("data/output/signals/$(filesuffix)/mid") && mkpath("data/output/signals/$(filesuffix)/mid")
!ispath("data/output/signals/$(filesuffix)/high") && mkpath("data/output/signals/$(filesuffix)/high")

## Read in data
rawdata = DataFrame(load("data/raw/sgs_data.dta"))
rawdata[:date] = makedate.(rawdata[:date])

## Set data range to use in samples
startindex = findfirst(isequal(Date(1980, 1)), rawdata[:date])
endindex = findfirst(isequal(Date(2017, 12)), rawdata[:date])

## Estimate model for each horizion we want with no signals
results = sampleSignals(Vector{Float64}(rawdata[series]), Vector{Date}(rawdata[:date]), 1, 12:12, startindex:endindex;
    D = 3, burnin = 10_000, Nrun = 10_000, initialburn = 100_000, initialNrun = 1, signalLen = 0, noise = 0.5, 
    noiseSamples= 0)
saveresults(results, "data/output/$(filesuffix)/") 

## Estimate model for each horizon we want with signals for different noise levels
results = sampleSignals(Vector{Float64}(rawdata[series]), Vector{Date}(rawdata[:date]), 1, 12:12, startindex:12:endindex;
    D = 3, burnin = 10_000, Nrun = 10_000, initialburn = 100_000, initialNrun = 1, signalLen = 1, noise = 0.5, 
    noiseSamples= 300)
saveresults(results, "data/output/signals/$(filesuffix)/low/") 

results = sampleSignals(Vector{Float64}(rawdata[series]), Vector{Date}(rawdata[:date]), 1, 12:12, startindex:12:endindex;
    D = 3, burnin = 10_000, Nrun = 10_000, initialburn = 100_000, initialNrun = 1, signalLen = 1, noise = 1.5, 
    noiseSamples= 300)
saveresults(results, "data/output/signals/$(filesuffix)/mid/") 

results = sampleSignals(Vector{Float64}(rawdata[series]), Vector{Date}(rawdata[:date]), 1, 12:12, startindex:12:endindex;
    D = 3, burnin = 10_000, Nrun = 10_000, initialburn = 100_000, initialNrun = 1, signalLen = 1, noise = 3.0, 
    noiseSamples = 300)
saveresults(results, "data/output/signals/$(filesuffix)/high/") 

stateresults = smoothStates(Vector{Float64}(rawdata[series]), rawdata[:date], 1:endindex; D = 3, burnin = 100_000, Nrun = 10_000)
savesmoothresults(stateresults, "data/output/$(filesuffix)/")
