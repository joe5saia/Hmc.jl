###############################################################################
#= alt_forecasts.jl
Code to calculate forecasts from postiors of data
=##############################################################################
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
using CSV, DataFrames, Statistics, LinearAlgebra


## Read in data
dfmeans = CSV.read(joinpath(root_dir, "data/output/official/filtered_means_summary.csv"))
dfstates = CSV.read(joinpath(root_dir, "data/output/official/filtered_state_probs_summary.csv"))
dfA = CSV.read(joinpath(root_dir, "data/output/official/filtered_trans_probs_summary.csv"))

dfforecasts = DataFrame(date = dfmeans[!,:date], forecast_post = fill(1.0, size(dfmeans,1)), forecast_proper = fill(1.0, size(dfmeans,1)) )

for i in 1:size(dfmeans,1)
    μ = collect(dfmeans[i,2:4])
    π = collect(dfstates[i,2:4])
    A = reshape(collect(dfA[i,2:end]), 3,3)
    dfforecasts[i,:forecast_post] = dot(π' * A^12, μ)
end

dfforecastsorg = CSV.read(joinpath(root_dir, "data/output/official/forecasts_summary.csv"))
dfforecasts[!, :forecast_proper] = dfforecastsorg[!, :forecast_12_mean]
CSV.write(joinpath(root_dir, "data/output/official/forecast_comparision.csv"), dfforecasts)