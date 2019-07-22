##############################################################################
#= 
Program: dispersion.jl
Author: Joe Saia
Date: July 2019
Info: Calculates the mean, std, and IQR of the postior forecasts from the 
Monte Carlo simulations with signals
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

for s in ["official", "alter"]
    for n in ["1.0", "3.0", "7.0"]
        Hmc.calcdispersion(joinpath(root_dir, "data/output/signals/$(s)/noise_$(n)"))
    end
end
