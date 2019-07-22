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
if ispath("/moto/sscc/projects/biasedexpectations")
    root_dir = "/moto/sscc/projects/biasedexpectations"
elseif ispath("/research/hmc")
    root_dir = "/research/hmc"
else 
    @error "No valid directory for root directory found"
    exit(1)
end
cd(root_dir)


## Supply the data series as a command line argument 
if length(ARGS) !== 2
    @error "Must supply either \"official\" or \"alter\" as command line argument and date index as second argument e.g. \n\tjulia --project code/run_hmm.jl official 244"
    exit(1)
end

if ARGS[1] == "official"
    const dataseries = :offic_inf
    const series = "official"
    println("Using official inflation")
elseif ARGS[1] == "alter"
    const dataseries = :alter_inf
    const series = "alter"
    println("Using alternative inflation")
else 
    @error "Must supply either \"official\" or \"alter\" as first command line argument and date index as second argument e.g. \n\tjulia --project code/run_hmm.jl official 244"
    exit(2)
end

dateindex = parse(Int, ARGS[2])
println("End date index is $(dateindex)")


## Load packages 
using Pkg
Pkg.activate(root_dir)
push!(LOAD_PATH, joinpath(root_dir, "src"))
using Hmc

using StatFiles
using DataFrames
using Dates


## Read in data
rawdata = DataFrame(load("data/raw/sgs_data.dta"))
rawdata[:date] = makedate.(rawdata[:date])

## Set data range to use in samples
startindex = findfirst(isequal(Date(1970, 1)), rawdata[:date])

## Check to make sure that supplied end date index is valid
if length(rawdata[:date]) < dateindex
    @error "End date index is too high. Last date is $(rawdata[end,:date]) at index $(length(rawdata[:date])). You supplied $(dateindex)"
    exit(4)
elseif dateindex < startindex
    @error "End date index is too low. Start index is $(startindex). You supplied $(dateindex)"
    exit(5)
end

## Run parameters
Nrun = 5
burnin = 10_000
signalNrun = 5
signalburnin = 1_000
signalLen = 1
D = 3
noiseSamples = 3

## Estimate model for each horizion we want with no signals
println("Estimating base model without Signals. Date range is $(rawdata[startindex, :date]) to $(rawdata[dateindex, :date])")
!ispath("data/output/$(series)/") && mkpath("data/output/$(series)/")
estimateSignal(Vector{Float64}(rawdata[dataseries]), Vector{Date}(rawdata[:date]); startIndex=startindex, endIndex=dateindex, horizons=12,
    D = D, burnin = burnin, Nrun = Nrun, signalburnin = signalburnin, signalNrun = signalNrun, signalLen = 0, noise = 0.0, 
    noiseSamples = noiseSamples, σsignal = Float64[], savenosignal= true, series = series)


## Estimate model for each horizon we want with signals for different noise levels
for noiselevel in [1.0 3.0 7.0]
    println("Running noise = $(noiselevel) sample")
    !ispath("data/output/signals/$(series)/noise_$(noiselevel)") && mkpath("data/output/signals/$(series)/noise_$(noiselevel)")
    estimateSignal(Vector{Float64}(rawdata[dataseries]), Vector{Date}(rawdata[:date]); startIndex=startindex, endIndex=dateindex, horizons=12,
        D = D, burnin = burnin, Nrun = Nrun, signalburnin = signalburnin, signalNrun = signalNrun, signalLen = 1, noise = noiselevel, 
        noiseSamples = noiseSamples, σsignal= Float64[], savenosignal = false, series = series)
end
