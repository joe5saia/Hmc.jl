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
@everywhere root_dir = "/moto/sscc/projects/biasedexpectations/"
@everywhere cd(root_dir)

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

## Run parameters
Nrun = 3_000
burnin = 1_000
initialNrun = 1
initialburn = 100_000
signalLen = 1
D = 3
noiseSamples = 300


## Estimate model for each horizion we want with no signals
results = sampleSignals(Vector{Float64}(rawdata[series]), Vector{Date}(rawdata[:date]), 1, 12:12, startindex:endindex;
    D = D, burnin = burnin, Nrun = Nrun, initialburn = initialburn, initialNrun = initialNrun, signalLen = 0, noiseSamples = 3,
    noise = 0.0)
saveresults(results, "data/output/$(filesuffix)/") 

## Estimate model for each horizon we want with signals for different noise levels
println("Running low noise sample")
results = sampleSignals(Vector{Float64}(rawdata[series]), Vector{Date}(rawdata[:date]), 1, 12:12, startindex:endindex;
    D = D, burnin = burnin, Nrun = Nrun, initialburn = initialburn, initialNrun = initialNrun, signalLen = signalLen, noiseSamples = noiseSamples,
    noise = 1.0)
saveresults(results, "data/output/signals/$(filesuffix)/low/") 


println("Running mid noise sample")
results = sampleSignals(Vector{Float64}(rawdata[series]), Vector{Date}(rawdata[:date]), 1, 12:12, startindex:endindex;
    D = D, burnin = burnin, Nrun = Nrun, initialburn = initialburn, initialNrun = initialNrun, signalLen = signalLen, noiseSamples = noiseSamples,
    noise = 3.0)
saveresults(results, "data/output/signals/$(filesuffix)/mid/") 

println("Running high noise sample")
results = sampleSignals(Vector{Float64}(rawdata[series]), Vector{Date}(rawdata[:date]), 1, 12:12, startindex:endindex;
    D = D, burnin = burnin, Nrun = Nrun, initialburn = initialburn, initialNrun = initialNrun, signalLen = signalLen, noiseSamples = noiseSamples,
    noise = 7.0)
saveresults(results, "data/output/signals/$(filesuffix)/high/") 

println("Running smooth sample")
stateresults = smoothStates(Vector{Float64}(rawdata[series]), rawdata[:date], 1:endindex; D = 3, burnin = burnin, Nrun = Nrun)
savesmoothresults(stateresults, "data/output/$(filesuffix)/")
