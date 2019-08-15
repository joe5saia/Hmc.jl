if ispath("/moto/sscc/projects/biasedexpectations")
    root_dir = "/moto/sscc/projects/biasedexpectations"
elseif ispath("/research/hmc")
    root_dir = "/research/hmc"
else 
    @error "No valid directory for root directory found"
    exit(1)
end
cd(root_dir)
using StatFiles, DataFrames, Dates, CSV

using Pkg
Pkg.activate(root_dir)
push!(LOAD_PATH, joinpath(root_dir, "src"))
using Hmc

## Read in data
rawdata = DataFrame(load("data/raw/sgs_data.dta"))
rawdata[!, :date] = makedate.(rawdata[!, :date])
CSV.write("data/processed/inflation.csv", rawdata)