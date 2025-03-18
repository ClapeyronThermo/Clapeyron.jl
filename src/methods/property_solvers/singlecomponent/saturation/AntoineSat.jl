"""
    AntoineSaturation <: SaturationMethod
    AntoineSaturation(;T0 = nothing,
                        vl = nothing,
                        vv = nothing,
                        f_limit = 0.0,
                        atol = 1e-8,
                        rtol = 1e-12,
                        max_iters = 10^4,
                        crit = nothing,
                        crit_retry = false)

Saturation method for `saturation_temperature` .Default method for saturation temperature from Clapeyron 0.3.7. It solves the Volume-Temperature system of equations for the saturation condition.

If only `T0` is provided, `vl` and `vv` are obtained via [`x0_sat_pure`](@ref). If `T0` is not provided, it will be obtained via [`x0_saturation_temperature`](@ref). It is recommended to overload `x0_saturation_temperature`, as the default starting point calls [`crit_pure`](@ref), resulting in slower than ideal times.
`f_limit`, `atol`, `rtol`, `max_iters` are passed to the non linear system solver.
"""
struct AntoineSaturation{T,V,C} <: SaturationMethod
    T0::T
    vl::V
    vv::V
    f_limit::Float64
    atol::Float64
    rtol::Float64
    max_iters::Int
    crit::C
    crit_retry::Bool
end

function AntoineSaturation(;T0 = nothing,
    vl = nothing,
    vv = nothing,
    f_limit = 0.0,
    atol = 1e-8,
    rtol = 1e-12,
    max_iters = 10^2,
    crit = nothing,
    crit_retry = false)
    C = typeof(crit)
    if T0 === vl === vv === nothing
        AntoineSaturation{Nothing,Nothing,C}(nothing,nothing,nothing,f_limit,atol,rtol,max_iters,crit,crit_retry)
    elseif !(T0 === nothing) & (vl === vv === nothing)
        return AntoineSaturation{typeof(T0),Nothing,C}(T0,vl,vv,f_limit,atol,rtol,max_iters,crit,crit_retry)
    elseif (T0 === nothing) & !(vl === nothing) & !(vv === nothing)
        vl,vv = promote(vl,vv)
        return AntoineSaturation{Nothing,typeof(vl),C}(T0,vl,vv,f_limit,atol,rtol,max_iters,crit,crit_retry)
    elseif !(T0 === nothing) & !(vl === nothing) & !(vv === nothing)
        T0,vl,vv = promote(T0,vl,vv)
        return AntoineSaturation{typeof(T0),typeof(vl),C}(T0,vl,vv,f_limit,atol,rtol,max_iters,crit,crit_retry)
    else
        throw(error("invalid specification of AntoineSaturation"))
    end
end

function Obj_Sat_Temp(model::EoSModel, T, V_l, V_v,p,scales)
    fun(_V) = eos(model, _V, T,SA[1.])
    A_l,Av_l = Solvers.f∂f(fun,V_l)
    A_v,Av_v =Solvers.f∂f(fun,V_v)
    g_l = muladd(-V_l,Av_l,A_l)
    g_v = muladd(-V_v,Av_v,A_v)
    (p_scale,μ_scale) = scales
    F1 = -(Av_l+p)*p_scale
    F2 = -(Av_v+p)*p_scale
    F3 = (g_l-g_v)*μ_scale
    return SVector(F1,F2,F3)
end

function saturation_temperature_impl(model,p,method::AntoineSaturation{TT,VV,CC}) where {TT,VV,CC}
    R̄ = Rgas(model)
    ps,μs = equilibria_scale(model)
    if isnothing(method.T0)
        T0,Vl,Vv = x0_saturation_temperature(model,p)
        if !(isnothing(method.vl) && isnothing(method.vv))
            Vl,Vv = method.vl,method.vv
        end
    elseif isnothing(method.vl) && isnothing(method.vv)
        Vl,_ = x0_sat_pure(model,method.T0)
        T0 = method.T0
        Vv = Rgas(model)*T0/p
    else
        T0,Vl,Vv = method.T0,method.vl,method.vv
    end
    T0,Vl,Vv = promote(T0,Vl,Vv,oneunit(eltype(model)))
    nan = zero(T0)/zero(T0)
    fail = (nan,nan,nan)

    if isnan(T0)
        return fail
    end
    res,converged = try_2ph_pure_temperature(model,p,T0,Vl,Vv,ps,μs,method)
    converged && return res
    T2,_,_ = res

    if !isnothing(method.crit)
        Tc,pc,_ = method.crit
        p > pc && return fail
        T2 >= Tc && return fail
    end

    #try if we are too low or too high
    #one (or two) saturation pressure calculations are normally faster than a crit pure calculation
    (p2,vl2,vv2) = saturation_pressure(model,T2,ChemPotVSaturation(crit_retry = false))
    if !isnan(p2) #nice, psat(T2) exists, we can now produce a really good estimate of the saturation temperature
        dpdT = dpdT_saturation(model,vl2,vv2,T2)
        dTinvdlnp = -p2/(dpdT*T2*T2)
        Δlnp = log(p/p2)
        Tinv0 = 1/T2
        Tinv = Tinv0 + dTinvdlnp*Δlnp
        T3 = 1/Tinv
        (_,vl3,vv3) = saturation_pressure(model,T3,ChemPotVSaturation(crit_retry = false))
        res,converged = try_2ph_pure_temperature(model,p,T3,vl3,vv3,ps,μs,method)
        converged && return res
    end
    #no luck here, we need to calculate the critical point
    if !method.crit_retry
        return fail
    end

    if isnothing(crit)
        crit = crit_pure(model)
    end
    Tc,pc,vc = crit
    p > pc && return fail
    if 0.9 < p/pc < 1.0
        #you could still perform another iteration from a better initial point
        T3,Vl3,Vv3 = x0_saturation_temperature(model,p,crit)
        if !(Vl ≈ Vl3) && !(Vv ≈ Vv3) && !(T3 ≈ T0) #check if the initial points are not the same
            res,converged = try_2ph_pure_temperature(model,p,T3,Vl3,Vv3,ps,μs,method)
            converged && return res
        end
    end
    return fail
end

export AntoineSaturation
