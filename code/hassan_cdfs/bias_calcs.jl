##############################################################################
#= 
Program: calc_cdfs.jl
Author: Joe Saia
Date: November 2019
Info: Calculate a function for Hassan. See the outline directory for information
Output: Writes several csvs with postior statistics
=#
##############################################################################
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
using Hmc, Dates, CSV, Distributions, LinearAlgebra, Gadfly, DataFrames, Cairo, Statistics


#########################################################################
# TODO
# 1. Get forecast bias from filtered_forecats_xxxx.csv
# 2. Get 12-month ahead state probabilties from filtered_forecats_xxxx.csv
#    and means from filtered_means_xxxx.csv 
# 3. Do covariance calculation
#########################################################################

# Set up
dates = Date(1980,1):Dates.Year(1):Date(2018,1)
datadir = joinpath(root_dir, "data/output/official")

function calcbias(d)
  means = Matrix(CSV.read(joinpath(datadir, "filtered_means_$(d).csv"))[:, 2:end])
  totalbias = mean(Vector(CSV.read(joinpath(datadir, "forecasts_$(d).csv"))[:, 3]))
  states = Matrix(CSV.read(joinpath(datadir, "filtered_state_probs_$(d).csv"))[:, 2:end])
  trans = Matrix(CSV.read(joinpath(datadir, "filtered_trans_probs_$(d).csv"))[:,2:end])
  futstate = Array{Float64,2}(undef,size(states))
  for i in 1:size(futstate,1)
    futstate[i, :] = states[i, :]' * reshape(trans[i, :],3,3)^12
  end
  covv = cov(hcat(futstate, means))
  nab = mean(hcat(means, futstate), dims=1)
  uncbias = nab * covv * nab'
  return [uncbias[1], totalbias]
end

output = Array{Any,2}(undef, length(dates), 3)
for (i,d) in enumerate(dates)
  println(d)
  output[i,1] = d
  output[i,2:end] = calcbias(d)
end

df  = DataFrame(output, [:date, :uncerbias, :totalbias])
plot(stack(df,[:uncerbias, :totalbias]), x = :date, y = :value, color = :variable, Geom.line) |> PDF("code/hassan_cdfs/outline/biasestime.pdf")