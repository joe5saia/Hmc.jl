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

for s in ["official"]
    Hmc.runaggregate(joinpath(root_dir, "data/output/$(s)"))
    for n in ["0.01"], m in ["12"]
        if ispath(joinpath(root_dir, "data/output/signals_$(s)_noise_$(n)_len_$(m)"))
            Hmc.runaggregate(joinpath(root_dir, "data/output/signals_$(s)_noise_$(n)_len_$(m)"))
        else
            println("$(joinpath(root_dir, "data/output/signals_$(s)_noise_$(n)_len_$(m)")) is not a valid path" )
        end
    end
end




