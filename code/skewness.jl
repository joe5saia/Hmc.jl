###############################################################################
#= skewness.jl
Code to calculate skweness of final period inflation distribution for Laura 
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
using CSV, DataFrames, Statistics, Distributions

## Read in data
dfmeans = CSV.read(joinpath(root_dir, "data/output/official/filtered_means_summary.csv"))
dfvars = CSV.read(joinpath(root_dir, "data/output/official/filtered_variances_summary.csv"))
dftrans = CSV.read(joinpath(root_dir, "data/output/official/filtered_trans_probs_summary.csv"))
dfstates = CSV.read(joinpath(root_dir, "data/output/official/filtered_state_probs_summary.csv"))

## Set the date to calculate skewness for
ldate = maximum(dfmeans[!,:date])
println("using data from $ldate")

## Make distribution that describes the mixture distribution inflation follows
f = MixtureModel(Normal[
    Normal(dfmeans[dfmeans[!,:date] .== ldate, :state_1_mean][1], sqrt(dfvars[dfvars[!,:date] .== ldate, :state_1_mean][1])),
    Normal(dfmeans[dfmeans[!,:date] .== ldate, :state_2_mean][1], sqrt(dfvars[dfvars[!,:date] .== ldate, :state_2_mean][1])),
    Normal(dfmeans[dfmeans[!,:date] .== ldate, :state_3_mean][1], sqrt(dfvars[dfvars[!,:date] .== ldate, :state_3_mean][1])) ],
    collect(dfstates[dfstates[!,:date] .== ldate, [:state_1_mean,:state_2_mean, :state_3_mean] ][1,:]) )

## draw from mixture model and then caluclate sample skewness
data = rand(f,10_000_000)
skew = sum(((data .- mean(f))./std(f)).^3)/length(data)
println("Skewness is $skew")