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
elseif ispath("/app")
    root_dir = "/app"
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
        sampleRange=startindex:dateindex, 
        signalRange=2:1, 
        endIndex = dateindex,
        horizons=[12],
        D=3, 
        burnin = 10, 
        Nrun = 10, 
        signalburnin = 3, 
        signalNrun = 5, 
        noiseSamples = 4,
        series = series
    )
else
    ## Run parameters
    println("Using actual settings\n")
    opt = Hmc.estopt(
        Vector{Float64}(rawdata[!,dataseries]), 
        Vector{Date}(rawdata[!,:date]),
        sampleRange=startindex:dateindex, 
        signalRange=2:1, 
        endIndex=dateindex,
        horizons=[12],
        D = 3,
        burnin = 100_000,
        Nrun = 250_000,
        signalburnin = 100_000,
        signalNrun = 250_000,
        noiseSamples = 100,
        series = series
    )
end


if false
    #Estimate model for each horizion we want with no signals
    println("Estimating base model without Signals. Date range is $(Hmc.startdate(opt)) to $(Hmc.enddate(opt))")
    p = "data/output/$(opt.series)/"
    !ispath(p) && mkpath(p)
    Random.seed!(opt.seed)
    samples = Hmc.estimatemodel(opt)
    Hmc.saveresults(samples, opt, p; hassignals = false)
end
if false
    #Estimate model for each horizon we want with signals for different noise levels
    signalLen = 1
    opt.sampleRange = startindex:dateindex+signalLen
    opt.signalRange = dateindex+1:dateindex+signalLen
    opt.signalSave = dateindex+1:dateindex+signalLen
    Hmc.update_itators!(opt)
    for noiselevel in [0.1 1.0 3.0]
        println("Running noise = $(noiselevel) sample")
        opt.noise=noiselevel
        opt.σsignal=0
        p = "data/output/signals_$(opt.series)_noise_$(opt.noise)_len_$(signalLen)"
        !ispath(p) && mkpath(p)
        Random.seed!(opt.seed)
        samples = Hmc.estimatesignals!(opt)
        Hmc.saveresults(samples, opt, p; hassignals = true)
    end
end
if false
    signalLen = 12
    opt.sampleRange = startindex:dateindex+signalLen
    opt.signalRange = dateindex+1:dateindex+signalLen
    opt.signalSave = dateindex+1:dateindex+signalLen
    Hmc.update_itators!(opt)
    for noiselevel in [0.1 1.0 3.0]
        println("Running noise = $(noiselevel) sample")
        opt.noise=noiselevel
        opt.σsignal=0
        p = "data/output/signals_$(opt.series)_noise_$(opt.noise)_len_$(signalLen)"
        !ispath(p) && mkpath(p)
        Random.seed!(opt.seed)
        samples = Hmc.estimatesignals!(opt)
        Hmc.saveresults(samples, opt, p; hassignals = true)
    end
end

if true
    ## Make everything a signal 
    opt.sampleRange = startindex:dateindex
    opt.signalRange = startindex:dateindex
    opt.signalSave = dateindex-1:dateindex
    Hmc.update_itators!(opt)
    #for noiselevel in [0.1 0.3 1.0]
    for noiselevel in [0.6]
        println("Running noise = $(noiselevel) sample")
        opt.noise=noiselevel
        opt.σsignal=0
        p = "data/output/signals_$(opt.series)_noise_$(opt.noise)_allsignal"
        !ispath(p) && mkpath(p)
        Random.seed!(opt.seed)
        samples = Hmc.estimatesignals!(opt)
        Hmc.saveresults(samples, opt, p; hassignals = true)
    end

end
