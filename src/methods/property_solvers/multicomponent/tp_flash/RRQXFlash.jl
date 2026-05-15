#hooks to PTFlashWrapper
update_pressure!(model,p) = nothing
update_temperature!(model,T) = nothing

function update_K_QX!(model,p,T,w,β,phases,non_inw,vec_cache,dlnϕ_cache,verbose_cache)
    x,y,z =  w
    K,lnK = vec_cache
    #using cache
    v1 = lnK[1]
    v2 = lnK[2]
    phasex,phasey = phases
    non_inx,non_iny = non_inw

    lnϕx, volx = modified_lnϕ(model, p, T, x, dlnϕ_cache; phase = phasex)
    lnK .= lnϕx
    lnϕy, voly = modified_lnϕ(model, p, T, y, dlnϕ_cache; phase = phasey)

    lnK .-= lnϕy
    K .= exp.(lnK)
    lnK[1] = volx
    lnK[2] = voly
    for i in 1:length(K)
        non_inx[i] && (K[i] = Inf)
        non_iny[i] && (K[i] = 0)
    end
    x,y = update_rr!(K,β,z,x,y,non_inx,non_iny,false)

    OF = sum(x) - sum(y)
    verbose,iter_count,is_temperature = verbose_cache
    if verbose     
        iter_count[] += 1
        iter = iter_count[]
        if is_temperature
            X = T
        else
            X = p
        end
        @info "$(__pad_val(iter,4))    $(__pad_val(X,16))   $(__pad_val(OF,16)) $(repr(K,context = :compact => true))"
    end
    return OF

end

function update_K_QT!(logp,params)
    p = exp(logp)
    model,β,T,w,phases,non_inw,vec_cache,dlnϕ_cache,verbose = params
    return update_K_QX!(model,p,T,w,β,phases,non_inw,vec_cache,dlnϕ_cache,verbose)
end

function update_K_QP!(Tinv,params)
    T = 1/Tinv
    model,β,p,w,phases,non_inw,vec_cache,dlnϕ_cache,verbose = params
    update_temperature!(model,T)
    return update_K_QX!(model,p,T,w,β,phases,non_inw,vec_cache,dlnϕ_cache,verbose)
end

"""
    RRQXFlash{T}(;kwargs...)

Method to solve non-reactive multicomponent, two-phase flash problem, using a generalized formulation.

Only two phases are supported. If `K0` is `nothing`, it will be calculated via fugacity coefficients at p,T conditions.

### Keyword Arguments:
- `equilibrium` (optional) = equilibrium type ":vle" for liquid vapor equilibria, ":lle" for liquid liquid equilibria, `:unknown` if not specified
- `p0` (optional), initial guess pressure, ignored if pressure is one of the flash specifications.
- `T0` (optional), initial guess temperature, ignored if temperature is one of the flash specifications.
- `K0` (optional), initial guess for the K-values.
- `x0` (optional), initial guess for the composition of phase x.
- `y0` = optional, initial guess for the composition of phase y.
- `vol0` = optional, initial guesses for phase x and phase y volumes.
- `atol` = absolute tolerance to stop the calculation.
- `rtol` = relative tolerance to stop the calculation.
- `max_iters` = maximum number of iterations
- `flash_result::FlashResult`: can be provided instead of `x0`,`y0` and `vol0` for initial guesses.
"""
struct RRQXFlash{P,T} <: FlashMethod
    equilibrium::Symbol
    T0::Union{P,Nothing}
    p0::Union{P,Nothing}
    K0::Union{Vector{T},Nothing}
    x0::Union{Vector{T},Nothing}
    y0::Union{Vector{T},Nothing}
    v0::Union{Tuple{T,T},Nothing}
    atol::Float64
    rtol::Float64
    max_iters::Int
    verbose::Bool
end

function Solvers.primalval(method::RRQXFlash{P,T}) where {P,T}
    if P == Nothing
        λP = Nothing
    else
        λP = Solvers.primal_eltype(P)
    end
    
    if T == Nothing
        λT = Nothing
    else
        λT = Solvers.primal_eltype(P)
    end
    return RRQXFlash{λP,λT}(method.equilibrium,primalval(method.T0),primalval(method.p0),primalval(method.K0),primalval(method.x0),primalval(method.y0),primalval(method.v0),method.atol,method.rtol,method.max_iters,method.verbose)
end

Base.eltype(method::RRQXFlash{T}) where T = T

function index_reduction(m::RRQXFlash,idx::AbstractVector)
    equilibrium,T0,p0,K0,x0,y0,v0,atol,rtol,max_iters = m.equilibrium,m.T0,m.p0,m.K0,m.x0,m.y0,m.v0,m.atol,m.rtol,m.max_iters
    K0 !== nothing && (K0 = K0[idx])
    x0 !== nothing && (x0 = x0[idx])
    y0 !== nothing && (y0 = y0[idx])
    return RRQXFlash(;equilibrium,T0,p0,K0,x0,y0,v0,atol,rtol,max_iters,method.verbose)
end

index_reduction(m::RRQXFlash{Nothing,Nothing},idx::AbstractVector) = m

numphases(::RRQXFlash) = 2

function RRQXFlash(;equilibrium = :unknown,
                        T0 = nothing,
                        p0 = nothing,
                        K0 = nothing,
                        x0 = nothing,
                        y0 = nothing,
                        v0 = nothing,
                        rtol = 1e-14,
                        atol = 1e-12,
                        max_iters = 100,
                        flash_result = nothing,
                        verbose = false,)
    !(is_vle(equilibrium) | is_lle(equilibrium) | is_unknown(equilibrium))  && throw(error("invalid equilibrium specification for RRQXFlash"))
    if flash_result isa FlashResult
        comps,β,volumes = flash_result.compositions,flash_result.fractions,flash_result.volumes
        np = numphases(flash_result)
        np != 2 && incorrect_np_flash_error(RRQXFlash,flash_result)
        w1,w2 = comps[1],comps[2]
        v = (volumes[1],volumes[2])
        P00 = flash_result.data.p
        T00 = flash_result.data.T
        return RRQXFlash(;equilibrium = equilibrium,T0 = T00,p0 = P00,x0 = w1,y0 = w2,v0 = v,rtol = rtol,atol = atol,max_iters = max_iters,verbose = verbose)
    end

    if K0 == x0 == y0 === nothing #nothing specified
        #is_lle(equilibrium)
        T = Nothing
    else
        if !isnothing(K0) & isnothing(x0) & isnothing(y0) #K0 specified
            T = eltype(K0)
        elseif isnothing(K0) & !isnothing(x0) & !isnothing(y0)  #x0, y0 specified
            T = eltype(x0)
        else
            throw(error("invalid specification of initial points"))
        end
    end

    if T == Nothing && v0 !== nothing
        TT = Base.promote_eltype(v0[1],v0[2])
        _v0 = (v0[1],v0[2])
    elseif T != nothing && v0 !== nothing
        TT = Base.promote_eltype(one(T),v0[1],v0[2])
        _v0 = (v0[1],v0[2])
    else
        TT = T
        _v0 = v0
    end

    if T0 === nothing && p0 === nothing
        S = Nothing
    elseif T0 !== nothing && p0 !== nothing
        S = typeof(T0*p0)
    else
        S = typeof(something(T0,p0))
    end
    return RRQXFlash{S,TT}(equilibrium,T0,p0,K0,x0,y0,_v0,atol,rtol,max_iters,verbose)
end

function qp_flash_impl(model,β,p,z,method::RRQXFlash)
    cached_model = __tpflash_cache_model(model,p,NaN,z,method.equilibrium)
    if cached_model != model
        return qp_flash_impl(cached_model,β,p,z,method)
    end
    flash0 = qp_flash_x0(model,β,p,z,method)
    TT = Base.promote_eltype(model,β,p,z)
    x = flash0.compositions[1]
    y = flash0.compositions[2]
    K = similar(z,TT)
    lnK = similar(K)
    x .= flash0.compositions[1]
    y .= flash0.compositions[2]
    w = (x,y,z)
    if is_vle(method.equilibrium) || is_unknown(method.equilibrium)
        phasex,phasey = :liquid,:vapour
    elseif is_lle(method.equilibrium)
        phasex,phasey = :liquid,:liquid
    end
    phases = (phasex,phasey)
    nc = length(model)
    #TODO: support this
    non_inx = FillArrays.Fill(false,nc)
    non_iny = FillArrays.Fill(false,nc)
    non_inw = (non_inx,non_iny)

    dlnϕ_cache = ∂lnϕ_cache(model, β, p, x, Val{false}())
    vec_cache = (K,lnK)
    #model,p,T,w,β,Kx,phases,non_inw,vec_cache,dlnϕ_cache,spec
    verbose = get_verbosity(method)
    verbose_cache = (get_verbosity(method),Ref(0),true)
    params = (model,β,p,w,phases,non_inw,vec_cache,dlnϕ_cache,verbose_cache)
    Tinv0 = 1/temperature(flash0)
    prob = Roots.ZeroProblem(update_K_QP!,Tinv0)
    
    verbose && @info "_______________________________________________________
      iter    T                  OF               K"

    Tinv = Roots.solve(prob,Roots.Order1(),params,atol = method.atol,rtol = method.rtol)
    T = 1/Tinv

@info "________________________________________________________________________________________
      Final K values:        $K
      Final temperature:     $T

"
    n = sum(z)
    resize!(lnK,2)
    resize!(K,2)
    K[1] = 1 - β
    K[2] = β
    x ./= sum(x)
    y ./= sum(y)
    K .*= n
    lnK ./= n
    lle = is_liquid(phasex) && is_liquid(phasey)
    vapour_idx = lle ? -1 : 2
    return FlashResult(flash0.compositions,K,lnK,FlashData(p,T,0.0,vapour_idx))
end

function qt_flash_impl(model::M,β,T,z,method::RRQXFlash) where M
    cached_model = __tpflash_cache_model(model,β,T,z,method.equilibrium)
    if cached_model != model
        return qt_flash_impl(cached_model,β,T,z,method)
    end
    flash0 = qt_flash_x0(model,β,T,z,method)
    TT = Base.promote_eltype(model,β,T,z)
    x = flash0.compositions[1]
    y = flash0.compositions[2]
    K = similar(z,TT)
    lnK = similar(K)
    x .= flash0.compositions[1]
    y .= flash0.compositions[2]
    w = (x,y,z)
    if is_vle(method.equilibrium) || is_unknown(method.equilibrium)
        phasex,phasey = :liquid,:vapour
    elseif is_lle(method.equilibrium)
        phasex,phasey = :liquid,:liquid
    end
    phases = (phasex,phasey)
    nc = length(model)
    #TODO: support this
    non_inx = FillArrays.Fill(false,nc)
    non_iny = FillArrays.Fill(false,nc)
    non_inw = (non_inx,non_iny)

    dlnϕ_cache = ∂lnϕ_cache(model, β, T, x, Val{false}())
    vec_cache = (K,lnK)
    #model,p,T,w,β,Kx,phases,non_inw,vec_cache,dlnϕ_cache,spec
    verbose = get_verbosity(method)
    verbose_cache = (verbose,Ref(0),false)
    params = (model,β,T,w,phases,non_inw,vec_cache,dlnϕ_cache,verbose_cache)
    logp0 = log(pressure(flash0))
    prob = Roots.ZeroProblem(update_K_QT!,logp0)

    verbose && @info "_______________________________________________________
      iter    p                  OF               K"
    logp = Roots.solve(prob,Roots.Order1(),params,atol = method.atol,rtol = method.rtol)
    p = exp(logp)

@info "________________________________________________________________________________________
      Final K values:        $K
      Final pressure:        $p

"
    n = sum(z)
    resize!(lnK,2)
    resize!(K,2)
    K[1] = 1 - β
    K[2] = β
    x ./= sum(x)
    y ./= sum(y)
    K .*= n
    lnK ./= n

    if model isa PTFlashWrapper
        if iszero(lnK[1]) && is_liquid(phasex)
            lnK[1] = volume(model,p,T,x,phase = phasex)
        end

        if iszero(lnK[2]) && is_liquid(phasey)
            lnK[2] = volume(model,p,T,y,phase = phasey)
        end
    end
    lle = is_liquid(phasex) && is_liquid(phasey)
    vapour_idx = lle ? -1 : 2
    return FlashResult(flash0.compositions,K,lnK,FlashData(p,T,0.0,vapour_idx))
end

function bubble_pressure_impl(model::EoSModel,T,z,method::RRQXFlash)
    return qp_to_bubblep(model,T,z,method)
end

function dew_pressure_impl(model::EoSModel,T,z,method::RRQXFlash)
    return qp_to_dewp(model,T,z,method)
end

function bubble_temperature_impl(model::EoSModel,p,z,method::RRQXFlash)
    return qp_to_bubblep(model,p,z,method)
end

function dew_temperature_impl(model::EoSModel,p,z,method::RRQXFlash)
    return qp_to_dewp(model,p,z,method)
end

export RRQXFlash