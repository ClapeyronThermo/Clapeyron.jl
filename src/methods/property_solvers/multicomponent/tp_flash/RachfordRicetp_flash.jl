"""
    RRTPFlash{T}(;kwargs...)

Method to solve non-reactive multicomponent flash problem by Rachford-Rice equation.

Only two phases are supported. if `K0` is `nothing`, it will be calculated via the Wilson correlation.

### Keyword Arguments:
- `equilibrium`: `:vle` for liquid vapor equilibria, `:lle` for liquid liquid equilibria, `:unknown` if not specified
- `K0`: initial guess for the K-values.
- `x0`: initial guess for the composition of phase x.
- `y0`: initial guess for the composition of phase y.
- `vol0`: initial guesses for phase x and phase y volumes.
- `K_tol`: tolerance to stop the calculation.
- `ss_iters`: number of Successive Substitution iterations to perform.
- `nacc`: accelerate successive substitution method every nacc steps. Should be a integer bigger than 3. Set to 0 for no acceleration.
- `second_order`: wheter to solve the gibbs energy minimization using the analytical hessian or not.
- `noncondensables`: arrays with names (strings) of components non allowed on the liquid phase. In the case of LLE equilibria, corresponds to the `x` phase.
- `nonvolatiles`: arrays with names (strings) of components non allowed on the vapour phase. In the case of LLE equilibria, corresponds to the `y` phase.
- `flash_result::FlashResult`: can be provided instead of `x0`,`y0` and `vol0` for initial guesses.

"""
struct RRTPFlash{T} <: TPFlashMethod
    equilibrium::Symbol
    K0::Union{Vector{T},Nothing}
    x0::Union{Vector{T},Nothing}
    y0::Union{Vector{T},Nothing}
    v0::Union{Tuple{T,T},Nothing}
    K_tol::Float64
    max_iters::Int
    nacc::Int
    noncondensables::Union{Nothing,Vector{String}}
    nonvolatiles::Union{Nothing,Vector{String}}
end

Base.eltype(method::RRTPFlash{T}) where T = T

function index_reduction(m::RRTPFlash,idx::AbstractVector)
    equilibrium,K0,x0,y0,v0,K_tol,max_iters,nacc,noncondensables,nonvolatiles = m.equilibrium,m.K0,m.x0,m.y0,m.v0,m.K_tol,m.max_iters,m.nacc,m.noncondensables,m.nonvolatiles
    K0 !== nothing && (K0 = K0[idx])
    x0 !== nothing && (x0 = x0[idx])
    y0 !== nothing && (y0 = y0[idx])
    return RRTPFlash{eltype(m)}(equilibrium,K0,x0,y0,v0,K_tol,max_iters,nacc,noncondensables,nonvolatiles)
end

function RRTPFlash(;equilibrium = :unknown,
    K0 = nothing,
    x0 = nothing,
    y0 = nothing,
    v0 = nothing,
    K_tol = 1e-10,
    max_iters = 100,
    nacc = 5,
    noncondensables = nothing,
    nonvolatiles = nothing,
    flash_result = nothing)
    ss_iters = max_iters
    if flash_result isa FlashResult
        np = numphases(flash_result)
        np != 2 && incorrect_np_flash_error(RRTPFlash,flash_result)
    end
    
    nonvolatiles isa String && (nonvolatiles = [nonvolatiles])
    noncondensables isa String && (noncondensables = [noncondensables])

    #we call Michelsen to check if the arguments are correct.
    m = MichelsenTPFlash(;equilibrium,K0,x0,y0,v0,K_tol,ss_iters,nacc,noncondensables,nonvolatiles,flash_result)
    return RRTPFlash{eltype(m)}(m.equilibrium,m.K0,m.x0,m.y0,m.v0,m.K_tol,max_iters,m.nacc,m.noncondensables,m.nonvolatiles)
end

function tp_flash_impl(model::EoSModel, p, T, z, method::RRTPFlash)

    model_cached = __tpflash_cache_model(model,p,T,z,method.equilibrium)

    x,y,β,v = tp_flash_michelsen(model_cached,p,T,z;equilibrium = method.equilibrium,
    K0 = method.K0, x0 = method.x0, y0 = method.y0, vol0 = method.v0,
    K_tol = method.K_tol,itss = method.max_iters, nacc=method.nacc,
    non_inx_list=method.noncondensables, non_iny_list=method.nonvolatiles,
    reduced = true, use_opt_solver = false)

    g = __tpflash_gibbs_reduced(model_cached,p,T,x,y,β,method.equilibrium)
    comps = [x,y]
    volumes = [v[1],v[2]]
    βi = [1-β ,β]
    return FlashResult(comps,βi,volumes,FlashData(p,T,g))
end

export RRTPFlash
