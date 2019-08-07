using Dates, Test, Statistics, Random, BenchmarkTools

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

#Number of observations
T = 500
# State Transition Matrix
A99 = [0.50 0.50; 0.20 0.80]
# means and variances
μ99 = [-5.0, 4.0]
σ99 = [1.0, 0.50]

Y99, X99 = generateData(A99, μ99, σ99, T, 2)
dates = range(Date(1970, 1, 1), step = Dates.Month(1), length = T)

opt = Hmc.estopt(
       Vector{Float64}(Y99), 
       Vector{Date}(dates),
       sampleRange=1:T-24,
       signalRange=2:1,
       endIndex=T-24,
       horizons=[12],
       D = 2,
       burnin = 3_000,
       Nrun = 1_000,
       signalburnin = 1,
       signalNrun = 1,
       series = "test"
)



Random.seed!(opt.seed)
samples = Hmc.estimatemodel(opt)


println("True mean: $(μ99)\nEstimated mean $(vec(mean(samples.μ, dims=1)))")
println("True variances: $(σ99)\nEstimated variances $(vec(mean(samples.σ, dims=1)))")
println("True A: $(vec(A99))\nEstimated A $(vec(mean(samples.A, dims=1)))")

@test all(isapprox.(vec(mean(samples.μ, dims=1)), μ99, atol=0.3))
@test all(isapprox.(vec(mean(samples.σ, dims=1)), σ99, atol=0.5))

exit()

Random.seed!(opt.seed)
samples = Hmc.estimatesignals!(opt)
println("True mean: $(μ99)\nEstimated mean $(vec(mean(samples.μ, dims=1)))")
println("True variances: $(σ99)\nEstimated variances $(vec(mean(samples.σ, dims=1)))")

@test all(isapprox.(vec(mean(samples.μ, dims=1)), μ99, atol=0.3))
@test all(isapprox.(vec(mean(samples.σ, dims=1)), σ99, atol=0.5))
