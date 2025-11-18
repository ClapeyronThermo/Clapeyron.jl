## Bubble pressure solver
"""
    ActivityBubblePressure(kwargs...)

Function to compute [`bubble_pressure`](@ref) using Activity Coefficients.
On activity coefficient models it solves the problem via succesive substitucion.
On helmholtz-based models, it approximates the activity coefficient using the saturated pure state as reference.

Inputs:
- `y0 = nothing`: optional, initial guess for the vapor phase composition
- `p0 = nothing`: optional, initial guess for the bubble pressure `[Pa]`
- `vol0 = nothing`: optional, initial guesses for the liquid and vapor phase volumes `[m³]`
- `atol = 1e-8`: optional, absolute tolerance of the non linear system of equations
- `rtol = 1e-12`: optional, relative tolerance of the non linear system of equations
- `itmax_ss = 40`: optional, maximum number of sucesive substitution iterations
"""
struct ActivityBubblePressure{T} <: BubblePointMethod
    vol0::Union{Nothing,Tuple{T,T}}
    p0::Union{Nothing,T}
    y0::Union{Nothing,Vector{T}}
    nonvolatiles::Union{Nothing,Vector{String}}
    itmax_ss::Int64
    rtol_ss::Float64
end

function ActivityBubblePressure(;vol0 = nothing,
                                p0 = nothing,
                                y0 = nothing,
                                nonvolatiles = nothing,
                                itmax_ss = 40,
                                rtol_ss = 1e-8)

    if p0 == y0 == vol0 == nothing
        return ActivityBubblePressure{Nothing}(vol0,p0,y0,nonvolatiles,itmax_ss,rtol_ss)
    elseif (p0 == y0 == nothing) && !isnothing(vol0)
        vl,vv = promote(vol0[1],vol0[2])
        return ActivityBubblePressure{typeof(vl)}(vol0,p0,y0,nonvolatiles,itmax_ss,rtol_ss)
    elseif (vol0 == y0 == nothing) && !isnothing(p0)
        p0 = float(p0)
        return ActivityBubblePressure{typeof(p0)}(vol0,p0,y0,nonvolatiles,itmax_ss,rtol_ss)
    elseif (p0 == vol0 == nothing) && !isnothing(y0)
        T = eltype(y0)
        return ActivityBubblePressure{T}(vol0,p0,y0,nonvolatiles,itmax_ss,rtol_ss)
    elseif !isnothing(vol0) && !isnothing(p0) && !isnothing(y0)
        vl,vv,p0,_ = promote(vol0[1],vol0[2],p0,first(y0))
        T = eltype(vl)
        y0 = convert(Vector{T},y0)
        return ActivityBubblePressure{T}(vol0,p0,y0,nonvolatiles,itmax_ss,rtol_ss)
    elseif !isnothing(vol0) && !isnothing(y0)
        vl,vv,_ = promote(vol0[1],vol0[2],first(y0))
        T = eltype(vl)
        y0 = convert(Vector{T},y0)
        return ActivityBubblePressure{T}(vol0,p0,y0,nonvolatiles,itmax_ss,rtol_ss)
    elseif  !isnothing(p0) && !isnothing(y0)
        p0,_ = promote(p0,first(y0))
        T = eltype(p0)
        y0 = convert(Vector{T},y0)
        return ActivityBubblePressure{T}(vol0,p0,y0,nonvolatiles,itmax_ss,rtol_ss)
    else
        throw(error("invalid specification for bubble pressure"))
    end
end

function bubble_pressure_impl(model,T,x,method::ActivityBubblePressure)
    wrapper = PTFlashWrapper(model,NaN,T,x,:vle)
    return bubble_pressure_impl(wrapper,T,x,method)
end

function bubble_pressure_impl(model::CompositeModel,T,x,method::ActivityBubblePressure)
    return bubble_pressure_impl(model.fluid,T,x,method)
end

function bubble_pressure_impl(wrapper::PTFlashWrapper,T,x,method::ActivityBubblePressure)
    TT = Base.promote_eltype(wrapper,T,x)
    volatiles = comps_in_equilibria(component_list(model),method.nonvolatiles)
    p,vl,vv,y = bubble_pressure_init(model,T,x,method.vol0,method.p0,method.y0,volatiles)
    cache = ∂lnϕ_cache(model, p, T, x, Val{false}())
    lnK = similar(x,TT)
    pold = -p
    lnK = similar(x,TT)
    K = similar(x,TT)
    for k in 1:method.itmax_ss
        lnϕx, volx = modified_lnϕ(model, p, T, x, cache; phase = phasex, vol0 = volx)
        lnK .= lnϕx
        lnϕy, voly = modified_lnϕ(model, p, T, y, cache; phase = phasey, vol0 = voly)
        lnK .-= lnϕy
        lnK .+= log(p)
        K .= exp.(lnK)
        y .= K .* x
        pold = p
        p = sum(y)
        y ./= p
        err = abs(pold-p)/p
        !isfinite(err) && break
        err < method.rtol_ss && break
    end
    if iszero(vl)
        vl = volume(wrapper,p,T,x,phase = :liquid)
    end
    return (p,vl,vv,y)
end


export ActivityBubblePressure