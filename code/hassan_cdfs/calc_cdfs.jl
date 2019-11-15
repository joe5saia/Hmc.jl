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
using Hmc, Dates, CSV, Distributions, LinearAlgebra, Gadfly, DataFrames, Cairo

## Read in data
dates = Date(1980,1):Dates.Year(5):Date(2015,1)
datadir = joinpath(root_dir, "data/output/official")
ys = -5:.25:15
means = Matrix(CSV.read(joinpath(datadir, "filtered_means_$(dates[1]).csv"))[:, 2:end])
expectations = Array{Float64,2}(undef, size(means,1), length(ys))

for date in dates
  means = Matrix(CSV.read(joinpath(datadir, "filtered_means_$(date).csv"))[:, 2:end])
  vars = Matrix(CSV.read(joinpath(datadir, "filtered_variances_$(date).csv"))[:, 2:end])
  πs = Matrix(CSV.read(joinpath(datadir, "filtered_state_probs_$(date).csv"))[:, 2:end])
  for (j,y) = enumerate(ys), i = 1:size(means,1)
    expectations[i,j] = dot(cdf.(Normal(), (y .- means[i,:])./(sqrt.(vars[i,:]))), πs[i,:]) 
  end 
  expectationsbar = mean(expectations, dims=1)
  plotdata = DataFrame(x = ys, finverse = quantile.(Normal(), vec(expectationsbar)))
  ## Plot 
  plot(plotdata, layer(x=:finverse, y=:x, Geom.line),  Guide.title("Enddate = $(date)")) |> PDF(joinpath(root_dir, "code/hassan_cdfs", "$(date).pdf"))
end 