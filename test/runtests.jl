using Dates
using Test
include("../src/hmc.jl")
using .Hmc

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

results = Hmc.sampleAndForecastAll(Y99, dates, 1:T - 20, 0:2, T - 20:T - 20;
                    D = 2, burnin = 1_000, Nrun = 10_000,
                    initialburn = 10_000, initialNrun = 1)

@test all(isapprox.(results[2][2:end], μ99, atol=0.1))
@test all(isapprox.(results[3][2:end], σ99, atol=0.1))
