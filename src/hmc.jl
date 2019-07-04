module Hmc
export makedate, sampleAndForecastAll, smoothStates, saveresults, savesmoothresults, generateData

using Distributions
using Random
using LinearAlgebra
using StaticArrays
using StatFiles
using DataFrames
using Dates
using CSV

struct HyperParams{T}
    """
    Struct to hold hyper parameters/ model specifications that don't change
    """
    D::Int # Number of states
    N::Int # Number of observations
    M::Int # Number of observations that are signals
    κ::Float64 ## relative precision of signals 
    ξ::SVector{T,Float64} # Prior for means
    α::SVector{T,Float64} # sample size for prior of σ²  
    g::SVector{T,Float64} # sample size for β hyper prior
    h::SVector{T,Float64} # scale for for β hyper prior
    ν::SVector{T,Float64} # sample size for prior of means 
    function HyperParams{T}(D, N, M, κ, ξ, α, g, h, ν) where {T}
        new(D, N, M, κ, ξ, α, g, h, ν)
    end
end

"""
    HyperParams(Y::AbstractVector{Float64}, D::Int64)
    
Returns a HyperParam struct 
"""
function HyperParams(Y::AbstractVector{Float64}, D::Int64; M::Int64=0, κ::Float64=1)

    N = size(Y, 1)
    R = maximum(Y) - minimum(Y)
    #ξ = SVector{D}(range(median(Y) - 0.25*R, median(Y) + 0.25*R, length=D)) # small α means this doesn't matter much
    ξ = SVector{D}([mean(Y) for i in 1:D]) # small α means this doesn't matter much
    α = SVector{D}([1.0 for i in 1:D])
    g = SVector{D}([0.2 for i in 1:D])
    h = SVector{D}([10/R^2 for i in 1:D])
    ν = SVector{D}([0.1 for i in 1:D])
    return HyperParams{D}(D, N, M, κ, ξ, α, g, h, ν)
end

function makeParams(Y::AbstractArray{Float64,1}, D::Int; μ₀=Float64[], A₀=Float64[], σ₀=Float64[], ρ₀=Float64[], β₀=Float64[])
    """
    Function to create new parameters and return them
    Initialize parameters to something reasonable to minimize problmatic edge cases in burn in
    Means are spaced evenly over data range, 
    variances equal to square root of overall data variance 
    Probability objects choosen to be uniform 
    States are selected so that the observation is closest to state mean 
    """
    N = size(Y,1)
    A = defaults(A₀, "A", ones(D,D)./D, D^2)
    for i in 1:D
        !isprobvec(A[i,:]) && error("$i row of A is not a probability")
    end
    R = maximum(Y) - minimum(Y)
    μ = defaults(μ₀, "μ", collect(range(median(Y) - 0.25*R, median(Y) + 0.25*R, length=D)) , D)
    σ = defaults(σ₀, "σ", repeat([std(Y)],D), D)
    ρ = defaults(ρ₀, "ρ", ones(D)./D, D)
    β = defaults(β₀, "β", repeat([1.0], D), D)
    π = ones(N,D)
    π2 = ones(N,D)
    P = ones(N,D,D)
    P2 = ones(N,D,D)
    X = ones(Int64,N)
    for i in 1:N
        X[i] = findmax(pdf.(Normal.(μ, σ), Y[i] ))[2]
    end
    for i in 1:N
        π[i,:] .= 1.0 / D
        π2[i,:] .= 1.0 / D
        P[i,:,:] .= 1.0 / D^2
        P2[i,:,:] .= 1.0 / D^2
    end
    return(A, ρ, μ, σ, β, π, π2, P, P2, X)
end

function defaults(X, varname, defaultValue, D)
    """
    Helper function to make sure everything is the right size
    """
    if length(X) == 0
        return defaultValue
    elseif length(X) != D
        error("$(varname) is length $(length(X)) but should be length $(D)")
    else
        return X
    end
end

function generateData(A, μ, σ, T::Int, D::Int)
    """
    Generates test data. Return the observable and states
    """
    Random.seed!(123)     
    Q = [Distributions.Categorical(A[k,:]) for k in 1:D]
    X = Array{Int64,1}(undef, T)
    X[1] = 1
    for t in 2:T
        X[t] = rand(Q[X[t - 1]])
    end

    # Simulate Observables
    Φ = vec(mapslices(z->Distributions.Normal(z[1], z[2]), [μ sqrt.(σ)], dims = 2))
    Y = Array{Float64,1}(undef, T)
    for t in 1:T
        Y[t] = rand(Φ[X[t]])
    end
    return Y,X
end

function update_μσ!(μ::AbstractVector{Float64}, σ::AbstractVector{Float64}, β::AbstractVector{Float64}, X::AbstractVector{Int64}, Y::AbstractVector{Float64}, hp::HyperParams)
    """
    Updates mean and variance in place 
    Draw σ^2 from inverseΓ(a,b) then draw  μ|σ from N(m,s)
    """
    N = hp.N
    M = hp.M
    κ = hp.κ
    D = hp.D
    ν = hp.ν
    ξ = hp.ξ
    α = hp.α
    Ni = zeros(eltype(X),D) # Number of states in i
    Mi = zeros(eltype(X),D) # number of signals in i
    S = zeros(eltype(Y),D) # sum of observables in i
    Sm = zeros(eltype(Y),D) # sum of signals in i
    S2 = zeros(eltype(Y),D) # sum of (observables - ybar)^2 in i
    Sm2 = zeros(eltype(Y),D) # sum of (signals - sbar)^2 in i
    ## Statistics 
    for t in 1:N-M
        i = X[t]
        Ni[i] += oneunit(i)
        S[i] += Y[t]
    end
    ybar = S./Ni
    if M > 0
        for t in 1+N-M:N
            i = X[t]
            Mi[i] += oneunit(i)
            Sm[i] += Y[t]
        end
    end
    sbar = Sm./Mi
    totalbar = @. (S + Sm)/(Ni + Mi)
    for t in 1:N-M
        i = X[t]
        S2[i] += (Y[t] - ybar)^2
    end
    if M > 0
        for t in 1+N-M:N
            i = X[t]
            Sm2[i] += (Y[t] - sbar)^2
        end
    end
    Meff = Mi./(1+κ) # effective number of signals
    Neff = Ni .+ Meff # effective number of total observations
    a = @. α + 0.5*Ni + 0.5*Mi
    b = @. β + 0.5*S2 + 0.5/(1+κ)*Sm2 + 0.5*Neff*ν/(Neff + ν)*(totalbar - ξ).^2
    for i in 1:D
        σ[i] = rand( Distributions.InverseGamma(a[i], b[i]) )
    end
        m = @. (S + sbar * Meff + ν*ξ)/(Neff + ν)
        s =  @. sqrt(σ/(Neff + ν))
    for i in 1:D
        μ[i] =  rand( Distributions.Normal(m[i], s[i]) )
    end
end

function update_β!(β::AbstractVector{Float64}, σ::AbstractVector{Float64}, hp::HyperParams)
    """
    Update the β hyper parameter in place
    """
    for i in 1:hp.D
        a = hp.g[i] + hp.α[i]
        b = hp.h[i] + σ[i]
        β[i] = rand(Distributions.InverseGamma(a, b))
    end
end

function update_ρ!(ρ::AbstractVector{Float64}, X::AbstractVector{Int64}, hp::HyperParams, Y::Vector{Float64})
    """
    Update the initial state parameter in place
    """
    a = ones(hp.D)
    ρ[:] = rand(Distributions.Dirichlet(a))
end

function update_A!(A::AbstractMatrix{Float64}, X::AbstractVector{Int64}, hp::HyperParams, Y::Vector{Float64})
    """
    Update the transistion probability matrix in place
    """
    Trans = ones(Int64,hp.D, hp.D)
    for i in 1:hp.N-1
       Trans[X[i],X[i+1]] += 1
    end
    for i in 1:hp.D
        A[i,:] = rand(Distributions.Dirichlet(Trans[i,:]))
    end
end

function forwardupdate_P!(P::AbstractArray{Float64,3}, π::AbstractMatrix{Float64}, A::AbstractMatrix{Float64}, μ::AbstractVector{Float64}, σ::AbstractVector{Float64}, ρ::AbstractVector{Float64}, hp::HyperParams, Y::Vector{Float64})
    """
    Forward recurision for the state probabilities in place
    """
    π[:] .= 0.0 
    dists = Array{Normal{Float64}}(undef,hp.D)
    for s in 1:hp.D
        dists[s] = Normal(μ[s],sqrt(σ[s]))
    end
    total = 0.0
    for r in 1:hp.D, s in 1:hp.D
        P[1,r,s] = ρ[r] * A[r,s] * pdf(dists[s], Y[1])
        total += P[1,r,s]
    end
    for r in 1:hp.D, s in 1:hp.D
        P[1,r,s] /= total
        π[1,s] += P[1,r,s]
    end

    for t in 2:hp.N
        total = 0.0
        for r in 1:hp.D, s in 1:hp.D
            P[t,r,s] = π[t-1,r] * A[r,s] * pdf(dists[s], Y[t])
            total += P[t,r,s]
        end
        for r in 1:hp.D, s in 1:hp.D
            P[t,r,s] /= total
            π[t,s] += P[t,r,s]
        end
        total = 0.0
        for s in 1:hp.D
            total += π[t,s]
        end
        for s in 1:hp.D
            π[t,s] /= total
        end
        isapprox(total, 0) && @warn "Forwardupdate_P has a near 0 total at step $(t). P is $(P[t-2,:,:])"
    end

end

function backwardupdate_P!(Pb::AbstractArray{Float64,3}, πb::AbstractMatrix{Float64}, Pf::AbstractArray{Float64,3}, πf::AbstractMatrix{Float64}, hp::HyperParams, Y::Vector{Float64})
    """
    Backward recurision for state probabilities in place
    """
    πb[:] .= 0.0
    Pb[end,:,:] = Pf[end,:,:]
    πb[end,:] = πf[end,:]
    for t in Iterators.reverse(1:hp.N-1)
        for r in 1:hp.D, s in 1:hp.D
            πb[t,r] += Pb[t+1,r,s]
        end
        for r in 1:hp.D, s in 1:hp.D
            Pb[t,r,s] = Pf[t,r,s] * πb[t,s]/πf[t,s]
        end
    end
end

function update_X!(X::AbstractVector{Int64}, πb::AbstractMatrix{Float64}, Pb::AbstractArray{Float64,3}, hp::HyperParams)
    """
    Update the state vector in place
    """
    p = zeros(hp.D)
    X[end] = rand(Distributions.Categorical(πb[end,:]))
    for k in Iterators.reverse(1:hp.N-1)
        total = 0.0
        for r in 1:hp.D
            p[r] = Pb[k+1,r,X[k+1]]
            total += p[r]
        end
        if total > eps()
            for j in 1:hp.D
                p[j] /= total
            end
        else
             for j in 1:hp.D
                p[j] = 1/hp.D
            end
        end
        X[k] = rand(Distributions.Categorical( p ))
    end
end

function gibbssweep!(μ::AbstractVector{Float64}, σ::AbstractVector{Float64}, β::AbstractVector{Float64}, ρ::AbstractVector{Float64}, A::AbstractMatrix{Float64}, πf::AbstractMatrix{Float64}, πb::AbstractMatrix{Float64}, Pf::AbstractArray{Float64,3}, Pb::AbstractArray{Float64,3}, X::AbstractVector{Int64}, hp::HyperParams, Y)
    """
    Run a Gibb's sweep and update all the Parameters
    After updating the P matrices, reorder everything so that the first state always has the lowest 
    mean, second state has the second highest and so forth. 
    After reordering, draw the new states
    All the parameters are updated in place
    """
    update_μσ!(μ, σ, β, X, Y, hp)
    update_β!(β, σ, hp::HyperParams)
    update_ρ!(ρ, X, hp::HyperParams, Y::Vector{Float64})
    update_A!(A, X, hp::HyperParams, Y::Vector{Float64})
    forwardupdate_P!(Pf, πf, A, μ, σ, ρ, hp::HyperParams, Y::Vector{Float64})
    backwardupdate_P!(Pb, πb, Pf, πf, hp::HyperParams, Y::Vector{Float64})
    # Order the states so that they are in increasing mean order
    order = sortperm(μ)
    μ[:] = μ[order]
    σ[:] = σ[order]
    β[:] = β[order]
    ρ[:] = ρ[order]
    Atmp = deepcopy(A)
    for (i1,i0) in enumerate(order)
        for (j1,j0) in enumerate(order)
            A[i1,j1] = Atmp[i0,j0] 
        end
    end
    πf[:] = πf[:, order]
    πb[:] = πf[:, order]
    update_X!(X, πb, Pb, hp)
end

function gibbssample!(A::AbstractMatrix{Float64}, 
                      ρ::AbstractVector{Float64},
                      μ::AbstractVector{Float64},
                      σ::AbstractVector{Float64},
                      β::AbstractVector{Float64},
                      πf::AbstractMatrix{Float64},
                      πb::AbstractMatrix{Float64},
                      Pf::AbstractArray{Float64,3},
                      Pb::AbstractArray{Float64,3},
                      X::AbstractVector{Int64},
                      hp::HyperParams,
                      Y::AbstractVector{Float64};
                      burnin::Int64 = 10, 
                      Nrun::Int64 = 2) 
    """
    Run the full sampling scheme and then calculate postior mean and return a named tuple with those parameters
    Takes the parameters as inputs and updates them in place
    """
    D = hp.D
    N = hp.N
    for i in 1:burnin
        gibbssweep!(μ, σ, β, ρ, A, πf, πb, Pf, Pb, X, hp::HyperParams, Y);
    end
    samples = Vector{Any}(undef,Nrun)
    for i in 1:Nrun 
        gibbssweep!(μ, σ, β, ρ, A, πf, πb, Pf, Pb, X, hp::HyperParams, Y);        
        samples[i] = (μ=μ, σ=σ, β=β, ρ=ρ, A=A, πf=πf, πb=πb)
    end
    ## Calculate the Postior 
    Apost = zeros(D,D)
    μpost = zeros(D)
    σpost = zeros(D)
    ρpost = zeros(D)
    πfpost = zeros(N,D)
    πbpost = zeros(N,D)
    for obs in samples
        Apost = Apost + obs.A
        μpost = μpost + obs.μ
        ρpost = ρpost + obs.ρ
        πfpost = πfpost + obs.πf
        πbpost = πbpost + obs.πb
        σpost = σpost + obs.σ
    end
    Apost = Apost / Nrun
    μpost = μpost / Nrun
    ρpost = ρpost / Nrun
    πfpost = πfpost / Nrun
    πbpost = πbpost / Nrun
    σpost = σpost / Nrun
    return (μ =μpost, σ = σpost, ρ = ρpost, A = Apost, πf = πfpost, πb= πbpost)
end

function makedate(x)
    """
    Helper function to make dates from data from stata
    Returns a Date
    """
    y = Int(x)
    year = y ÷ 12 + 1960
    month = mod(y, 12) + 1
    Dates.Date(year, month)
end

function sampleAndForecastAll(rawdata::AbstractVector{Float64}, dates::AbstractVector{Date}, dataRange, horizons, filterRange;
                     D::Int=2, burnin = 1_000, Nrun = 1_000, initialburn = 1_000, initialNrun = 1_000)
    """
    Function to run full HMM estimation for a range of samples
    Inputs -
    Rawdata: Vector of observables 
    dates: Vector of dates for each observable
    dataRange: Index range to select out observations to use in overall analysis
    horizons: Forecasts horizons to calculate
    filterRange: Index range to select the end point for each run. 
    ...
    Returns the forecasts, and postior means of the parameters for each run
    """
    ## Do a long initial burn in on first sample
    Y0 = rawdata[dataRange[1]:filterRange[1]]
    (A, ρ, μ, σ, β, πf, πb, Pf, Pb, X) = makeParams(Y0, D)
    hp0 = HyperParams(Y0, D)
    gibbssample!(A, ρ, μ, σ, β, πf, πb, Pf, Pb, X, hp0,
        Y0; burnin = initialburn, Nrun = initialNrun)
    
    μresults = Matrix{Any}(undef, length(filterRange),hp0.D+1)    
    σresults = Matrix{Any}(undef, length(filterRange),hp0.D+1)    
    Aresults = Matrix{Any}(undef, length(filterRange),hp0.D^2+1)    
    πbresults = Matrix{Any}(undef, length(filterRange),hp0.D+1)
    forecasts= Matrix{Any}(undef, length(filterRange), 2*length(horizons)+1)
    for (indx,j) in enumerate(filterRange)
        # Copy over parameters and states from previous filter and do short burn in
        println("Model run $j. End date: $(dates[j])")
        Y = rawdata[dataRange[1]:j]
        Xt = copy(X)
        (At, ρt, μt, σt, βt, πft, πbt, Pft, Pbt, X) = makeParams(Y, D)
        for i in 1:min(length(X), length(Xt)) # copy over states assuming same starting date
            X[i] = Xt[i]
        end
        hp = HyperParams(Y, D)
        parpost = gibbssample!(A, ρ, μ, σ, β, πft, πbt, Pft, Pbt, X, hp, Y;
                               burnin = burnin, Nrun = Nrun)
        ## Save results
        μresults[indx, 1] = dates[j]
        σresults[indx, 1] = dates[j]
        Aresults[indx, 1] = dates[j]
        πbresults[indx, 1] = dates[j]
        μresults[indx, 2:end] = parpost.μ
        σresults[indx, 2:end] = parpost.σ
        Aresults[indx, 2:end] = parpost.A[:]
        πbresults[indx, 2:end] = parpost.πb[end,:]
        forecasts[indx, 1] = dates[j]
        for (i,h) in enumerate(horizons)
            # println("Updating forecast horizon $h")
            forecasts[indx, 2*(i-1)+2] = forecast(parpost, hp0, h, rawdata[j+h])[1]
            forecasts[indx, 2*(i-1)+3] = forecast(parpost, hp0, h, rawdata[j+h])[2]
        end
    end
    return forecasts, μresults, σresults, Aresults, πbresults
end

function smoothStates(rawdata::AbstractVector{Float64}, dates::AbstractVector{Date}, dataRange;
                     D::Int=2, burnin = 1_000, Nrun = 1_000)
    """
    Function similiar to sampleAndForecastAll. Runs only one estimation and outputs the
    state probabilties for each day in the sample
    """
    Y0 = rawdata[dataRange]
    (A, ρ, μ, σ, β, πf, πb, Pf, Pb, X) = makeParams(Y0, D)
    hp0 = HyperParams(Y0, D)
    parpost = gibbssample!(A, ρ, μ, σ, β, πf, πb, Pf, Pb, X, hp0,
        Y0; burnin = burnin, Nrun = Nrun)
    
    πbresults = Matrix{Any}(undef, length(dataRange),hp0.D+1)
    for (indx,j) in enumerate(dataRange)
        πbresults[indx, 1] = dates[j]
        πbresults[indx, 2:end] = parpost.πb[j,:]
    end
    return πbresults
end

function forecast(par, hp::HyperParams, horizon, Yreal)
    """
    Calculate the expected mean horizon-periods ahead and the realized forecast error 
    """
    D = hp.D
    S0 = par.πb[end,:]'
    S1 = S0 * par.A^horizon 
    forecast = dot(S1,par.μ)
    forecasterror = forecast - Yreal
    return forecast, forecasterror
end

function saveresults(results, dir)
    """
    Write postiors to CSV's
    """
    forecasts, μresults, σresults, Aresults, πbresults = results

    forecastdf = DataFrame(forecasts)
    forecastnames = Vector{Symbol}(undef, size(forecasts,2))
    forecastnames[1] = :Date
    for j in 1:size(forecasts,2) ÷ 2
        forecastnames[2 * (j - 1) + 2] = Symbol("forecast_$j")
        forecastnames[2 * (j - 1) + 3] = Symbol("forecast_error_$j")
    end
    names!(forecastdf, forecastnames)
    Dname = [Symbol("state_$j") for j in 0:size(μresults,2)-1]
    Dname[1] = :Date
    D2name = [Symbol("transprob_$j") for j in 0:size(Aresults,2)-1]
    D2name[1] = :Date
    μresultsdf = DataFrame(μresults, Dname)
    σresultsdf = DataFrame(σresults, Dname)
    Aresultsdf = DataFrame(Aresults, D2name)
    πbresultsdf = DataFrame(πbresults, Dname)
    CSV.write(dir*"filtered_forecasts.csv", forecastdf)
    CSV.write(dir*"filtered_means.csv", μresultsdf)
    CSV.write(dir*"filtered_variance.csv", σresultsdf)
    CSV.write(dir*"filtered_trans_probs.csv", Aresultsdf)
    CSV.write(dir*"filtered_state_probs.csv", πbresultsdf)
end

function savesmoothresults(πbresults, dir)
    """
    Write the postior of the smoothed state probabilites to CSV
    """
    Dname = [Symbol("state_$j") for j in 0:size(πbresults,2)-1]
    Dname[1] = :Date
    πbresultsdf = DataFrame(πbresults, Dname)
    CSV.write(dir*"smoothed_state_probs.csv", πbresultsdf)
end

end # module
