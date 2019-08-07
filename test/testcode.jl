using Dates, Test, Statistics, Random, BenchmarkTools
using Distributions
using Random
using LinearAlgebra
using StaticArrays
using StatFiles
using DataFrames
using Dates
using CSV
using Distributed
using SharedArrays
using Glob


if ispath("/moto/sscc/projects/biasedexpectations")
    root_dir = "/moto/sscc/projects/biasedexpectations"
elseif ispath("/research/hmc")
    root_dir = "/research/hmc"
else 
    @error "No valid directory for root directory found"
    exit(1)
end
cd(root_dir)




mutable struct estopt
    rawdata::AbstractArray{Float64,1}
    dates::AbstractArray{Date,1}
    sampleRange::AbstractArray{Int,1}
    signalRange::AbstractArray{Int,1}
    signalSave::AbstractArray{Int,1}
    obsRange::AbstractArray{Int,1}
    sampleRangeit::Array{CartesianIndex{1},1}
    sampleRangeits2::Array{CartesianIndex{1},1}
    sampleRangeite2::Array{CartesianIndex{1},1}
    sampleRangeitback::Array{CartesianIndex{1},1}
    obsRangeit::Array{CartesianIndex{1},1}
    signalRangeit::Array{CartesianIndex{1},1}
    endIndex::Int
    horizons::AbstractArray{Int,1}
    D::Int
    burnin::Int
    Nrun::Int
    signalburnin::Int
    signalNrun::Int
    noise::Float64
    noiseSamples::Int
    σsignal::Float64
    series::String
    seed::Int
    function estopt(
        rawdata, 
        dates; 
        sampleRange=1:121, 
        signalRange=2:1, 
        signalSave=2:1,
        endIndex = 121,
        horizons=[12],
        D=3, 
        burnin = 1_000, 
        Nrun = 1_000, 
        signalburnin = 1_000, 
        signalNrun = 1_000, 
        noise = 0.0,
        noiseSamples = 1, 
        σsignal = 0.0,
        series = "offical",
        seed=1234
        )
        !issubset(signalRange, sampleRange) && @error "signalRange is not a subset of sampleRange"
        !issubset(signalSave, signalRange) && @error "signalSave is not a subset of signalRange"
        obsRange = setdiff(sampleRange, signalRange)
        sampleRangeit = [CartesianIndex(i) for i in sampleRange]
        sampleRangeits2 = [CartesianIndex(i) for i in sampleRange[2:end]] #range starting at 2
        sampleRangeite2 =  [CartesianIndex(i) for i in sampleRange[1:end-1]]#range ending 1 early
        sampleRangeitback = [CartesianIndex(i) for i in reverse(sampleRange[1:end-1])] #range starting from second to last going to 1
        obsRangeit = [CartesianIndex(i) for i in obsRange]
        signalRangeit = [CartesianIndex(i) for i in signalRange]
        new(rawdata, dates, sampleRange, signalRange, signalSave, obsRange, sampleRangeit, sampleRangeits2, sampleRangeite2, sampleRangeitback, obsRangeit, signalRangeit,
            endIndex, horizons, D, burnin, Nrun, signalburnin, signalNrun, noise, noiseSamples, σsignal, series, seed)
    end
end

function makey(x::estopt)
    return x.rawdata[x.sampleRange]
end

function makeysignals(x::estopt)
    return x.rawdata[x.signalRange]
end

function enddate(x::estopt, extra=0)
    return x.dates[x.endIndex+extra]
end

function startdate(x::estopt)
    return x.dates[first(x.sampleRange)]
end

function yobs(x::estopt, index)
    return x.rawdata[index]
end

function yend(x::estopt, extra=0)
    return x.rawdata[x.endIndex+extra]
end

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
function HyperParams(Y::AbstractVector{Float64}, D::Int64, M::Int64=0, κ::Float64=1.0)
    N = size(Y, 1)
    R = maximum(Y) - minimum(Y)
    #ξ = SVector{D}(range(median(Y) - 0.25*R, median(Y) + 0.25*R, length=D)) # small α means this doesn't matter much
    ξ = SVector{D}([mean(Y) for i in 1:D]) # small α means this doesn't matter much
    α = SVector{D}([1.0 for i in 1:D])
    g = SVector{D}([1.0 for i in 1:D])
    h = SVector{D}([10/R^(-2) for i in 1:D])
    ν = SVector{D}([1.0 for i in 1:D])
    return HyperParams{D}(D, N, M, κ, ξ, α, g, h, ν)
end

"""
Returns a HyperParam struct
"""
#function HyperParams(Y::AbstractVector{Float64}, D::Int64, M::Int64=0, κ::Float64=1.0)
function HyperParams(opt)
    Y = makey(opt)
    N = size(Y, 1)
    R = maximum(Y) - minimum(Y)
    #ξ = SVector{D}(range(median(Y) - 0.25*R, median(Y) + 0.25*R, length=D)) # small α means this doesn't matter much
    ξ = SVector{opt.D}([mean(Y) for i in 1:opt.D]) # small α means this doesn't matter much
    α = SVector{opt.D}([2.0 for i in 1:opt.D])
    g = SVector{opt.D}([2.0 for i in 1:opt.D])
    h = SVector{opt.D}([R^2 for i in 1:opt.D])
    ν = SVector{opt.D}([2.0 for i in 1:opt.D])
    return HyperParams{opt.D}(opt.D, N, length(opt.signalRange), opt.noise, ξ, α, g, h, ν)
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
        π[i,:] .= pdf.(Normal.(μ, σ), Y[i] )
        π2[i,:] = π[i,:]
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

function update_μσ!(μ::AbstractVector{Float64}, σ::AbstractVector{Float64}, β::AbstractVector{Float64}, X::AbstractVector{Int64}, Y::AbstractVector{Float64}, hp::HyperParams, opt::estopt)
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
    ybar = zeros(eltype(Y),D) # mean of observables in i
    sbar = zeros(eltype(Y),D) # mean of signals in i
    totalbar = zeros(eltype(Y),D) # mean of signals in i
    Sm = zeros(eltype(Y),D) # sum of signals in i
    S2 = zeros(eltype(Y),D) # sum of (observables - ybar)^2 in i
    Sm2 = zeros(eltype(Y),D) # sum of (signals - sbar)^2 in i
    ## Statistics
    ## Observations 
    for t in opt.obsRangeit
        i = X[t]
        Ni[i] += oneunit(i)
        S[i] += Y[t]
    end
    for i in 1:D
        if Ni[i] > 0
            ybar[i] = S[i]./Ni[i]
        else
            ybar[i] = zero(ybar[1])
        end
    end

    ## Signals
    for t in opt.signalRangeit
        i = X[t]
        Mi[i] += oneunit(i)
        Sm[i] += Y[t]
    end
    for i in 1:D
        if Mi[i] > 0
            sbar[i] = Sm[i]./Mi[i]
        else
            sbar[i] = zero(sbar[1])
        end
    end

    ## Averages
    for i in 1:D
        if (Ni[i] + Mi[i]) > 0
            totalbar[i] =  (S[i] + Sm[i])/(Ni[i] + Mi[i])
        else
            totalbar[i] = zero(sbar[1])
        end
    end

    ## Observations
    for t in opt.obsRangeit
        i = X[t]
        S2[i] += (Y[t] - ybar[i])^2
    end

    ## Signals 
    for t in opt.signalRangeit
        i = X[t]
        Sm2[i] += (Y[t] - sbar[i])^2
    end

    Meff = Mi./(1+κ) # effective number of signals
    Neff = Ni .+ Meff # effective number of total observations
    #println("\nNeff is $(Neff)")
    #println("Meff is $(Meff)")
    #println("\nNi is $(Ni)")
    #println("ybar is $(ybar)")
    #println("totalbar is $(totalbar)")
    #println("S2 is $(S2)")
    #println("Sm2 is $(Sm2)")
    #println("end term is $(@. 0.5*Neff*ν/(Neff + ν)*(totalbar - ξ)^2)")
    #println("Mi is $(Mi)")
    a = @. α + 0.5*Ni + 0.5*Mi
    b = @. β + 0.5*S2 + (0.5/(1+κ))*Sm2 + 0.5*Neff*ν/(Neff + ν)*(totalbar - ξ)^2
    #println("a is $(a)")
    #println("b is $(b)")
    #println("Expected Value is $(b./(a .- 1))")
    for i in 1:D
        try
            σ[i] = rand( Distributions.InverseGamma(a[i], b[i]) )
        catch
            println("a or b is not postive")
            println(b)
            println(β)
            println(S2)
            println(Sm2)
            println(Neff)
            println(totalbar)
        end
    end
        m = @. (S + Sm + ν*ξ)/(Neff + ν)
        s =  @. sqrt(σ/(Neff + ν))
    for i in 1:D
        μ[i] =  rand( Distributions.Normal(m[i], s[i]) )
    end
end

function update_β!(β::AbstractVector{Float64}, σ::AbstractVector{Float64}, hp::HyperParams)
    """
    Update the β hyper parameter in place
    """
    #for i in 1:hp.D
    #    a = hp.g[i] + hp.α[i]
    #    b = hp.h[i] + σ[i]^(-1)
    #    β[i] = rand(Distributions.Gamma(a, b))
    #end
    β .= 2.0
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

function forwardupdate_P!(P::AbstractArray{Float64,3}, π::AbstractMatrix{Float64}, A::AbstractMatrix{Float64}, μ::AbstractVector{Float64}, σ::AbstractVector{Float64}, ρ::AbstractVector{Float64}, hp::HyperParams, Y::Vector{Float64}, opt::estopt)
    """
    Forward recurision for the state probabilities in place
    """
    π[:] .= 0.0
    
    ## Create emission probability distributions
    distobs = Array{Normal{Float64}}(undef,hp.D)
    distsig = Array{Normal{Float64}}(undef,hp.D)
    for s in 1:hp.D
        distobs[s] = Normal(μ[s], sqrt(σ[s]))
        distsig[s] = Normal(μ[s], (1+hp.κ)*sqrt(σ[s]))
    end
    #println(distobs)

    total = 0.0
    if oneunit(Int) in opt.obsRange
        #println("in obs for 1")
        for s in 1:hp.D, r in 1:hp.D
            P[1,r,s] =  ρ[r] * A[r,s] * pdf(distobs[s], Y[1])
            total += P[1,r,s]
        end
    else
        #println("in sig for 1")
        for s in 1:hp.D, r in 1:hp.D
            P[1,r,s] = ρ[r] * A[r,s] * pdf(distsig[s], Y[1])
            total += P[1,r,s]
        end
    end
    for s in 1:hp.D, r in 1:hp.D
        P[1,r,s] /= total
        π[1,s] += P[1,r,s]
    end

    #println("pi 1 is $(π[1,:])")

    for (t,t2) in zip(opt.sampleRangeits2, opt.sampleRangeit)
        total = 0.0

        if t in opt.obsRangeit
            #println("in obs for $t")
            #println("pi is $(π[t2,:])")
            #println("Y is $(Y[t])")
            for s in 1:hp.D, r in 1:hp.D
                #println("Product is $(π[t2,r] * A[r,s] * pdf(distobs[s], Y[t]))")
                P[t,r,s] = π[t2,r] * A[r,s] * pdf(distobs[s], Y[t])
                #println("P is $(P[t,r,s])")
                total += P[t,r,s]
            end
            #println("P is $(P[t,:,:])")
        
        else
        
            #println("in sig for $t")
            for s in 1:hp.D, r in 1:hp.D
                P[t,r,s] = π[t2,r] * A[r,s] * pdf(distsig[s], Y[t])
                total += P[t,r,s]
           end

        end

        for s in 1:hp.D, r in 1:hp.D
            P[t,r,s] /= total
            π[t,s] += P[t,r,s]
        end

        #isapprox(total, 0) && @warn "Forwardupdate_P has a near 0 total at step $(t). P is $(P[t,:,:])"
    end
end

function backwardupdate_P!(Pb::AbstractArray{Float64,3}, πb::AbstractMatrix{Float64}, Pf::AbstractArray{Float64,3}, πf::AbstractMatrix{Float64}, hp::HyperParams, Y::Vector{Float64}, opt::estopt)
    """
    Backward recurision for state probabilities in place
    """
    πb[:,:] .= 0.0
    Pb[end,:,:] = Pf[end,:,:]
    πb[end,:] = πf[end,:]
    for t in opt.sampleRangeitback
        for s in 1:hp.D, r in 1:hp.D
            πb[t,r] += Pb[t[1]+1,r,s]
        end
        for s in 1:hp.D, r in 1:hp.D
            Pb[t,r,s] = Pf[t,r,s] * πb[t,s]/πf[t,s]
        end
    end
end

function update_X!(X::AbstractVector{Int64}, π::AbstractMatrix{Float64}, P::AbstractArray{Float64,3}, hp::HyperParams)
    """
    Update the state vector in place
    """
    #p = zeros(hp.D)
    #X[end] = rand(Distributions.Categorical(π[end,:]))
    for t in opt.sampleRangeitback
    #    total = 0.0
    #    for r in 1:hp.D
    #        p[r] = P[k+1,r,X[k+1]]
    #        total += p[r]
    #    end
    #    if total > eps()
    #        for j in 1:hp.D
    #            p[j] /= total
    #        end
    #    else
    #         for j in 1:hp.D
    #            p[j] = 1/hp.D
    #        end
    #    end
    #    X[k] = rand(Distributions.Categorical( p ))
    end
    return 1
end

function update_X!(X::AbstractVector{Int64}, π::AbstractMatrix{Float64}, P::AbstractArray{Float64,3}, hp::HyperParams, opt::estopt)
    """
    Update the state vector in place
    """    
    p = zeros(hp.D)
    X[end] = rand(Distributions.Categorical(π[end,:]))
    for k in opt.sampleRangeitback
        k2 = k[1]+1
        total = 0.0
        for r in 1:hp.D
            p[r] = P[k2,r,X[k2]]
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
    return 1
end



function gibbssweep!(μ::AbstractVector{Float64}, σ::AbstractVector{Float64}, β::AbstractVector{Float64}, ρ::AbstractVector{Float64}, A::AbstractMatrix{Float64}, πf::AbstractMatrix{Float64}, πb::AbstractMatrix{Float64}, Pf::AbstractArray{Float64,3}, Pb::AbstractArray{Float64,3}, X::AbstractVector{Int64}, hp::HyperParams, Y, opt::estopt)
    """
    Run a Gibb's sweep and update all the Parameters
    After updating the P matrices, reorder everything so that the first state always has the lowest
    mean, second state has the second highest and so forth.
    After reordering, draw the new states
    All the parameters are updated in place
    """
    update_μσ!(μ, σ, β, X, Y, hp, opt)
    update_β!(β, σ, hp::HyperParams)
    update_ρ!(ρ, X, hp::HyperParams, Y::Vector{Float64})
    update_A!(A, X, hp::HyperParams, Y::Vector{Float64})
    forwardupdate_P!(Pf, πf, A, μ, σ, ρ, hp::HyperParams, Y::Vector{Float64}, opt)
    backwardupdate_P!(Pb, πb, Pf, πf, hp::HyperParams, Y::Vector{Float64}, opt)
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
    πb[:] = πb[:, order]
    update_X!(X, πf, Pf, hp,opt)
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
                      Y::AbstractVector{Float64},
                      opt::estopt;
                      burnin::Int64 = 10,
                      Nrun::Int64 = 2
                      )
    """
    Run the full sampling scheme and then calculate postior mean and return a named tuple with those parameters
    Takes the parameters as inputs and updates them in place
    """
    D = hp.D
    N = hp.N
    ## Make some iterators that we will use alot

    #μ[:] = [-5.0, 4.0]
    #σ[:] = [1.0, 0.50]
    #ρ[:] = [0.5, 0.5]
    #A[:] = [0.50 0.50; 0.20 0.80]
    for i in 1:burnin
        gibbssweep!(μ, σ, β, ρ, A, πf, πb, Pf, Pb, X, hp::HyperParams, Y, opt);
    end

    ## do a gibbs sweep and save draws
    μsample = Array{Float64,2}(undef,Nrun, D)
    σsample = Array{Float64,2}(undef,Nrun, D)
    πbsample = Array{Float64,3}(undef,Nrun, N, D)
    Asample = Array{Float64,3}(undef,Nrun, D, D)
    for i in 1:Nrun
        gibbssweep!(μ, σ, β, ρ, A, πf, πb, Pf, Pb, X, hp::HyperParams, Y, opt);
        μsample[i, :] = μ
        σsample[i, :] = σ
        πbsample[i, :, :] = πb
        Asample[i, :, :] = A
    end
    return (μ = μsample, σ = σsample, πb = πbsample, A = Asample)
end


#Number of observations
T = 500
# State Transition Matrix
A99 = [0.50 0.50; 0.20 0.80]
# means and variances
μ99 = [-5.0, 4.0]
σ99 = [1.0, 0.50]

Y99, X99 = generateData(A99, μ99, σ99, T, 2)
dates = range(Date(1970, 1, 1), step = Dates.Month(1), length = T)

opt = estopt(
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



Y = makey(opt)
    (A, ρ, μ, σ, β, πf, πb, Pf, Pb, X) = makeParams(Y, opt.D)
    hp = HyperParams(Y, opt.D)


@benchmark update_μσ!(μ, σ, β, X, Y, hp, opt)
@benchmark update_β!(β, σ, hp::HyperParams)
@benchmark update_ρ!(ρ, X, hp::HyperParams, Y::Vector{Float64})
@benchmark update_A!(A, X, hp::HyperParams, Y::Vector{Float64})
@benchmark forwardupdate_P!(Pf, πf, A, μ, σ, ρ, hp::HyperParams, Y::Vector{Float64}, opt)
@benchmark backwardupdate_P!(Pb, πb, Pf, πf, hp::HyperParams, Y::Vector{Float64}, opt)
@benchmark update_X!(X, πf, Pf, hp)
@benchmark gibbssweep!(μ, σ, β, ρ, A, πf, πb, Pf, Pb, X, hp::HyperParams, Y, opt)


@code_warntype update_μσ!(μ, σ, β, X, Y, hp, opt)
@code_warntype update_β!(β, σ, hp::HyperParams)
@code_warntype update_ρ!(ρ, X, hp::HyperParams, Y::Vector{Float64})
@code_warntype update_A!(A, X, hp::HyperParams, Y::Vector{Float64})
@code_warntype forwardupdate_P!(Pf, πf, A, μ, σ, ρ, hp::HyperParams, Y::Vector{Float64}, opt)
@code_warntype backwardupdate_P!(Pb, πb, Pf, πf, hp::HyperParams, Y::Vector{Float64}, opt)



    ## Estimate model
 @benchmark gibbssample!(A, ρ, μ, σ, β, πf, πb, Pf, Pb, X, hp, Y, opt; burnin = opt.burnin, Nrun = opt.Nrun)
    forecasts= Array{Float64,2}(undef, opt.Nrun, 2*length(opt.horizons))