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
if length(ARGS) < 2
    @error "Must supply either \"official\" or \"alter\" as command line argument and date index as second argument e.g. \n\tjulia --project code/run_hmm.jl official 244"
    exit(1)
end

if ARGS[1] == "official"
    const dataseries = :offic_inf
    const series = "official"
    println("\n###############################################################################################\nUsing official inflation")
elseif ARGS[1] == "alter"
    const dataseries = :alter_inf
    const series = "alter"
    println("\n###############################################################################################\nUsing alternative inflation")
else 
    @error "Must supply either \"official\" or \"alter\" as first command line argument and date index as second argument e.g. \n\tjulia --project code/run_hmm.jl official 244"
    exit(2)
end

dateindex = parse(Int, ARGS[2])
println("End date index is $(dateindex)\n###############################################################################################")

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

## Check to make sure that supplied end date index is valid
if length(rawdata[!, :date]) < dateindex
    @error "End date index is too high. Last date is $(rawdata[end,:date]) at index $(length(rawdata[:date])). You supplied $(dateindex)"
    exit(4)
elseif dateindex < startindex
    @error "End date index is too low. Start index is $(startindex). You supplied $(dateindex)"
    exit(5)
end


if (length(ARGS) == 3) && (ARGS[3] == "test")
    ## Run parameters for testing
    println("Using testing run settings\n")
    
    opt = Hmc.estopt(
        Vector{Float64}(rawdata[!,dataseries]), 
        Vector{Date}(rawdata[!,:date]),
        startIndex=startindex,
        endIndex=dateindex,
        horizons=[12],
        D = 3,
        burnin = 300,
        Nrun = 10,
        signalburnin = 100,
        signalNrun = 3,
        signalLen = 0,
        noise = 0.0,
        noiseSamples = 3,
        σsignal = 0.0,
        savenosignal= true,
        series = series,
        seed = 1234,
    )
else
    ## Run parameters
    opt = Hmc.estopt(
        Vector{Float64}(rawdata[!,dataseries]), 
        Vector{Date}(rawdata[!,:date]),
        startIndex=startindex,
        endIndex=dateindex,
        horizons=[12],
        D = 3,
        burnin = 20_000,
        Nrun = 10_000,
        signalburnin = 10_000,
        signalNrun = 3_000,
        signalLen = 1,
        noise = 0.0,
        noiseSamples = 300,
        σsignal = 0.0,
        savenosignal= false,
        series = series,
        seed = 1234,
    )
end


# Estimate model for each horizion we want with no signals
println("Estimating base model without Signals. Date range is $(Hmc.startdate(opt)) to $(Hmc.enddate(opt))")
!ispath("data/output/$(opt.series)/") && mkpath("data/output/$(opt.series)/")
opt.signalLen=0
opt.savenosignal=true
Random.seed!(opt.seed)
samples = Hmc.estimatemodel(opt)
Hmc.saveresults(samples, opt; hassignals = false)



# Estimate model for each horizon we want with signals for different noise levels
opt.signalLen=1
opt.savenosignal=false
for noiselevel in [0.1 1.0 3.0]
    println("Running noise = $(noiselevel) sample")
    opt.noise=noiselevel
    p = "data/output/signals_$(opt.series)_noise_$(opt.noise)_len_$(opt.signalLen)"
    !ispath(p) && mkpath(p)
    Random.seed!(opt.seed)
    samples = Hmc.estimatesignals!(opt)
    Hmc.saveresults(samples, opt; hassignals = true)
end

opt.signalLen=1
opt.savenosignal=false
for noiselevel in [0.1 1.0 3.0]
    println("Running noise = $(noiselevel) sample")
    opt.noise=noiselevel
    p = "data/output/signals_$(opt.series)_noise_$(opt.noise)_len_$(opt.signalLen)"
    !ispath(p) && mkpath(p)
    Random.seed!(opt.seed)
    samples = Hmc.estimatesignals!(opt)
    Hmc.saveresults(samples, opt; hassignals = true)
end
