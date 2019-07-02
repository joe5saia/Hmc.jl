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

using StatFiles
using DataFrames
using Dates
include("../src/hmc.jl")
using .Hmc


rawdata = DataFrame(load("data/raw/sgs_data.dta"))
rawdata[:date] = makedate.(rawdata[:date])
startindex = findfirst(isequal(Date(1980, 01)), rawdata[:date])
endindex = findfirst(isequal(Date(2017, 12)), rawdata[:date])

results = sampleAndForecastAll(Vector{Float64}(rawdata[series]),
    Vector{Date}(rawdata[:date]), 1:endindex, 1:12, startindex:endindex;
    D = 3, burnin = 1_000, Nrun = 10_000, initialburn = 100_000, initialNrun = 1)
saveresults(results, "data/output/$(filesuffix)/") 

stateresults = smoothStates(Vector{Float64}(rawdata[series]), rawdata[:date], 1:endindex; D = 3, burnin = 10_000, Nrun = 10_000)
savesmoothresults(stateresults, "data/output/$(filesuffix)/")
