##############################################################################
#= 
Program: aggregate.jl
Author: Joe Saia
Date: July 2019
Info: Reads in output csvs from run_hmm.jl and calculates postior statistics
like the mean 
Output: Writes several csvs with postior statistics
=#
##############################################################################

using CSV, DataFrames, Dates, Statistics

## Supply the inflaton series as a command line argument 
if length(ARGS) == 0
    @error "Must supply either \"official\" or \"alter\" as command line argument e.g. \n\tjulia --project code/run_hmm.jl official"
    exit(1)
elseif ARGS[1] == "official"
    infseries = "official"
    println("Using official inflation")
elseif ARGS[1] == "alter"
    infseries = "alter"
    println("Using alternative inflation")
else 
    @error "Must supply either \"official\" or \"alter\" as command line argument e.g. \n\tjulia --project code/run_hmm.jl official"
    exit(2)
end

## Loop over files with gibbs draws we want to aggregate and do aggregation 
for noise in ["low", "mid", "high"]
    if ispath("data/output/signals/$(infseries)/$(noise)")
        for series in ["forecasts", "means", "state_probs", "trans_probs", "variance"]
            df = CSV.read("data/output/signals/$(infseries)/$(noise)/filtered_$(series).csv")
            df2 = aggregate(df, :Date, mean)
            CSV.write("data/output/signals/$(infseries)/$(noise)/filtered_$(series)_summary.csv", df2)
        end
    end
end

## Loop over files with gibbs draws we want to aggregate and do aggregation 
if ispath("data/output/$(infseries)")
    for series in ["forecasts", "means", "state_probs", "trans_probs", "variance"]
        df = CSV.read("data/output/$(infseries)/filtered_$(series).csv")
        df2 = aggregate(df, :Date, mean)
        CSV.write("data/output/$(infseries)/filtered_$(series)_summary.csv", df2)
    end
end







