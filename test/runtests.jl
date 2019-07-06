using Dates
using Test
using Statistics
using Distributed

push!(LOAD_PATH, "src")
using Hmc

#Number of observations
T = 5000
# State Transition Matrix
A99 = [0.25 0.75;
       0.60 0.40]
# means and variances
μ99 = [-5.0, 4.0]
σ99 = [1.0, 0.50]

Y99, X99 = generateData(A99, μ99, σ99, T, 2)
dates = range(Date(1970, 1, 1), step = Dates.Month(1), length = T)

results = sampleSignals(Y99, dates, 1, 8:12, T-20:T-20; D = 2, burnin = 10_000, Nrun = 10_000, 
                        initialburn = 10_000, initialNrun = 1, signalLen = 0, noise = 0.5, noiseSamples= 10)

forecasts, μresults, σresults, Aresults, πbresults, obsdates = results

@test all(isapprox.(vec(mean(μresults, dims=1)), μ99, atol=0.1))
@test all(isapprox.(vec(mean(σresults, dims=1)), σ99, atol=0.1))

println("True mean: $(μ99)\nEstimated mean $(vec(mean(μresults, dims=1)))")


results = sampleSignals(Y99, dates, 1, 8:12, T-20:T-20; D = 2, burnin = 1_000, Nrun = 1000, 
                        initialburn = 10_000, initialNrun = 1, signalLen = 1, noise = 0.5, noiseSamples= 10)

forecasts, μresults, σresults, Aresults, πbresults, obsdates = results
@test all(isapprox.(vec(mean(μresults, dims=1)), μ99, atol=0.1))
@test all(isapprox.(vec(mean(σresults, dims=1)), σ99, atol=0.1))
                        
