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

p = ARGS[1]
var = ARGS[2]

if ispath(p)
    Hmc.runaggregate(p, var)
else
    println(p * " is not a valid path" )
end
