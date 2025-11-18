"""
    ActivityDewPressure(kwargs...)

Function to compute [`dew_pressure`](@ref) using Activity Coefficients.
On activity coefficient models it solves the problem via succesive substitucion.
On helmholtz-based models, it uses the Chapman approximation for activity coefficients.

Inputs:
- `x0 = nothing`: optional, initial guess for the liquid phase composition
- `p0 = nothing`: optional, initial guess for the dew pressure `[Pa]`
- `vol0 = nothing`: optional, initial guesses for the liquid and vapor phase volumes `[m³]`
- `atol = 1e-8`: optional, absolute tolerance of the non linear system of equations
- `rtol = 1e-12`: optional, relative tolerance of the non linear system of equations
- `itmax_ss = 40`: optional, maximum number of sucesive substitution iterations
"""
struct ActivityDewPressure{T} <: DewPointMethod
    vol0::Union{Nothing,Tuple{T,T}}
    p0::Union{Nothing,T}
    x0::Union{Nothing,Vector{T}}
    noncondensables::Union{Nothing,Vector{String}}
    itmax_ss::Int64
    rtol_ss::Float64
end

function ActivityDewPressure(;vol0 = nothing,
                                p0 = nothing,
                                x0 = nothing,
                                noncondensables = nothing,
                                itmax_ss = 40,
                                rtol_ss = 1e-8)

    if p0 == x0 == vol0 == nothing
        return ActivityDewPressure{Nothing}(vol0,p0,x0,noncondensables,itmax_ss,rtol_ss)
    elseif (p0 == x0 == nothing) && !isnothing(vol0)
        vl,vv = promote(vol0[1],vol0[2])
        return ActivityDewPressure{typeof(vl)}(vol0,p0,x0,noncondensables,itmax_ss,rtol_ss)
    elseif (vol0 == x0 == nothing) && !isnothing(p0)
        p0 = float(p0)
        return ActivityDewPressure{typeof(p0)}(vol0,p0,x0,noncondensables,itmax_ss,rtol_ss)
    elseif (p0 == vol0 == nothing) && !isnothing(x0)
        T = eltype(x0)
        return ActivityDewPressure{T}(vol0,p0,x0,noncondensables,itmax_ss,rtol_ss)
    elseif !isnothing(vol0) && !isnothing(p0) && !isnothing(x0)
        vl,vv,p0,_ = promote(vol0[1],vol0[2],p0,first(x0))
        T = eltype(vl)
        x0 = convert(Vector{T},x0)
        return ActivityDewPressure{T}(vol0,p0,x0,noncondensables,itmax_ss,rtol_ss)
    elseif !isnothing(vol0) && !isnothing(x0)
        vl,vv,_ = promote(vol0[1],vol0[2],first(x0))
        T = eltype(vl)
        x0 = convert(Vector{T},x0)
        return ActivityDewPressure{T}(vol0,p0,x0,noncondensables,itmax_ss,rtol_ss)
    elseif  !isnothing(p0) && !isnothing(x0)
        p0,_ = promote(p0,first(x0))
        T = eltype(p0)
        x0 = convert(Vector{T},x0)
        return ActivityDewPressure{T}(vol0,p0,x0,noncondensables,itmax_ss,rtol_ss)
    else
        throw(error("invalid specification for dew pressure"))
    end
end

function dew_pressure_impl(model,T,y,method::ActivityDewPressure)
    wrapper = PTFlashWrapper(model,NaN,T,y,:vle)
    return dew_pressure_impl(wrapper,T,y,method)
end

function dew_pressure_impl(model::CompositeModel,T,y,method::ActivityDewPressure)
    return dew_pressure_impl(model.fluid,T,y,method)
end

function dew_pressure_impl(wrapper::PTFlashWrapper,T,y,method::ActivityDewPressure)
    condensables = comps_in_equilibria(component_list(model),method.noncondensables)
    p,volx,voly,y = bubble_pressure_init(model,T,x,method.vol0,method.p0,method.x0,condensables)
    cache = ∂lnϕ_cache(model, p, T, x, Val{false}())
    lnK = similar(x)
    pold = -p
    lnK = similar(x)
    K = similar(x)
    for k in 1:method.itmax_ss
        lnϕx, volx = modified_lnϕ(model, p, T, x, cache; phase = phasex, vol0 = volx)
        lnK .= lnϕx
        lnϕy, voly = modified_lnϕ(model, p, T, y, cache; phase = phasey, vol0 = voly)
        lnK .-= lnϕy
        lnK .+= log(p)
        K .= exp.(lnK)
        x .= y ./ K
        pold = p
        p = 1/sum(x)
        x .*= p
        err = abs(pold-p)/p
        !isfinite(err) && break
        err < method.rtol_ss && break
    end
    if iszero(vl)
        vl = volume(wrapper,p,T,x,phase = :liquid)
    end
    return (p,vl,vv,x)
end

export ActivityDewPressure