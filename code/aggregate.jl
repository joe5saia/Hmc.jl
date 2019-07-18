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
root_dir = "/moto/sscc/projects/biasedexpectations/"
#root_dir = "/research/hmc/"
cd(root_dir)

## Load packages 
using Pkg
Pkg.activate(root_dir)
push!(LOAD_PATH, "$(root_dir)src")
using Hmc

for s in ["offical", "alter"]
    aggregate(joinpath(root_dir, "data/output/$(s)"))
    for n in ["1.0", "3.0", "7.0"]
        aggregate(joinpath(root_dir, "data/output/signals/$(s)/noise_$n()"))
    end
end




