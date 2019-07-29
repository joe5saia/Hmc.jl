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

## Load packages 
using Pkg
Pkg.activate(root_dir)
push!(LOAD_PATH, joinpath(root_dir, "src"))
using Hmc

using StatFiles, DataFrames, Dates, Random


## Read in data
rawdata = DataFrame(load("data/raw/sgs_data.dta"))
rawdata[!, :date] = makedate.(rawdata[!, :date])

## Set data range to use in samples
startindex = findfirst(isequal(Date(1970, 1)), rawdata[!, :date])

    ## Run parameters
    opt = Hmc.estopt(
        Vector{Float64}(rawdata[!,:offic_inf]), 
        Vector{Date}(rawdata[!,:date]),
        startIndex=1,
        endIndex=576,
        horizons=[12],
        D = 3,
        burnin = 20_000,
        Nrun = 10_000,
        signalburnin = 10_000,
        signalNrun = 3_000,
        signalLen = 0,
        noise = 0.0,
        noiseSamples = 300,
        σsignal = 0.0,
        savenosignal= true,
        series = "official",
        seed = 1234,
    )


println("Estimating forecasts for each date. Date range is $(Hmc.startdate(opt)) to $(Hmc.enddate(opt))")
!ispath("data/output/$(opt.series)_insample/") && mkpath("data/output/$(opt.series)_insample/")
opt.savenosignal=true
Random.seed!(opt.seed)
samples = Hmc.estimatemodel(opt)
fcasts = Hmc.forecastinsample(samples, 12, opt)
Hmc.saveinsampleforecasts(fcasts, joinpath("data/output/$(opt.series)_insample", "forecats_insample.csv"))
