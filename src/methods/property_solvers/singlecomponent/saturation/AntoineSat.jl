"""
    AntoineSaturation <: SaturationMethod
    AntoineSaturation(;T0 = nothing,
                        vl = nothing,
                        vv = nothing,
                        f_limit = 0.0,
                        atol = 1e-8,
                        rtol = 1e-12,
                        max_iters = 10^4)
Saturation method for `saturation_temperature` .Default method for saturation temperature from Clapeyron 0.3.7. It solves the Volume-Temperature system of equations for the saturation condition.
    
If only `T0` is provided, `vl` and `vv` are obtained via [`x0_sat_pure`](@ref). If `T0` is not provided, it will be obtained via [`x0_saturation_temperature`](@ref). It is recommended to overload `x0_saturation_temperature`, as the default starting point calls [`crit_pure`](@ref), resulting in slower than ideal times.
`f_limit`, `atol`, `rtol`, `max_iters` are passed to the non linear system solver.
"""
struct AntoineSaturation{T,C} <: SaturationMethod
    T0::Union{Nothing,T}
    vl::Union{Nothing,T}
    vv::Union{Nothing,T}
    crit::C
    f_limit::Float64
    atol::Float64
    rtol::Float64
    max_iters::Int 
end

function AntoineSaturation(;T0 = nothing,
    vl = nothing,
    vv = nothing,
    crit = nothing,
    f_limit = 0.0,
    atol = 1e-8,
    rtol = 1e-12,
    max_iters = 10^2)
    C = typeof(crit)
    if T0 === vl === vv === nothing
        AntoineSaturation{Nothing,C}(nothing,nothing,nothing,crit,f_limit,atol,rtol,max_iters)
    elseif !(T0 === nothing) & (vl === vv === nothing)
        return AntoineSaturation{typeof(T0),C}(T0,vl,vv,crit,f_limit,atol,rtol,max_iters)
    elseif (T0 === nothing) & !(vl === nothing) & !(vv === nothing)
        vl,vv = promote(vl,vv)
        return AntoineSaturation{typeof(vl),C}(T0,vl,vv,crit,f_limit,atol,rtol,max_iters)
    elseif !(T0 === nothing) & !(vl === nothing) & !(vv === nothing)
        T0,vl,vv = promote(T0,vl,vv)
        return AntoineSaturation{typeof(vl),C}(T0,vl,vv,crit,f_limit,atol,rtol,max_iters)
    else
        throw(error("invalid specification of AntoineSaturation"))
    end
end

function Obj_Sat_Temp(model::EoSModel, F, T, V_l, V_v,p,scales,method::AntoineSaturation)
    fun(_V) = eos(model, _V, T,SA[1.])
    A_l,Av_l = Solvers.f∂f(fun,V_l)
    A_v,Av_v =Solvers.f∂f(fun,V_v)
    g_l = muladd(-V_l,Av_l,A_l)
    g_v = muladd(-V_v,Av_v,A_v)
    (p_scale,μ_scale) = scales
    F[1] = -(Av_l+p)*p_scale
    F[2] = -(Av_v+p)*p_scale
    F[3] = (g_l-g_v)*μ_scale
    return F
end

x0_saturation_temperature(model,p) = x0_saturation_temperature(model,p,AntoineSaturation())

function x0_saturation_temperature(model::EoSModel,p,::AntoineSaturation)
    coeffs = antoine_coef(model)
    coeffs === nothing && return x0_saturation_temperature(model,p,nothing)
    A,B,C = antoine_coef(model)
    lnp̄ = log(p / p_scale(model))
    T0 = T_scale(model)*(B/(A-lnp̄)-C)
    Vl,Vv = x0_sat_pure(model,T0)
    #take the liquid volume
    pl0 = pressure(model,Vl,T0)
    pv0 = pressure(model,Vv,T0)
    0.8*p < pl0 < 1.2*p && return (T0,Vl,Vv)
    #normally, the use of coefficients produce T0 values too far below from the real result, we add a correction
    R̄ = Rgas(model)
    Vv0 = R̄*T0/p
    ΔS = (VT_entropy(model,Vv0,T0) - VT_entropy(model,Vl,T0))
    dp = p - pl0
    Δv = Vv0 - Vl
    dpdt = ΔS/Δv
    dT = dp/dpdt    #(T - T0)
    T1 = T0 + dT
    if T1 < 0 #???
        T1 = T0 #error found while testing, better catch it earlier
    end
    Vv1 = 2*R̄*T1/p #add factor to allow easier iteration
    return (T1,Vl,Vv1)
end

#in case that there isn't any antoine coefficients:
#We aproximate to RK, use the cubic antoine, and perform refinement with one Clapeyron Saturation iteration 

function saturation_temperature_impl(model,p,method::AntoineSaturation)   
    R̄ = Rgas(model) 
    scales = scale_sat_pure(model)
    if isnothing(method.T0)
        T0,Vl,Vv = x0_saturation_temperature(model,p)
        if !(isnothing(method.vl) && isnothing(method.vv))
            Vl,Vv = method.vl,method.vv
        end
    elseif isnothing(method.vl) && isnothing(method.vv)
        Vl,Vv = x0_sat_pure(model,method.T0)
        T0 = method.T0
    else
        T0,Vl,Vv = method.T0,method.vl,method.vv
    end
    T0,Vl,Vv = promote(T0,Vl,Vv)
    nan = zero(T0)/zero(T0)
    fail = (nan,nan,nan)
    
    if isnan(T0)
        return fail
    end
    res,converged = try_sat_temp(model,p,T0,Vl,Vv,scales,method)
    converged && return res
    #it could be that the critical point isn't run there
    T2,_,_ = res
    
    crit = method.crit
    if !isnothing(crit)
        Tc,pc,_ = crit
        p > pc && return fail
        T2 >= Tc && return fail
    end

    #try if we are too low or too high
    #one (or two) saturation pressure calculations are normally faster than a crit pure calculation
    (p2,vl2,vv2) = saturation_pressure(model,T2,ChemPotVSaturation(crit_retry = false))
    if !isnan(p2) #nice, psat(T2) exists, we can now produce a really good estimate of the saturation temperature
        ΔHvap = (VT_enthalpy(model,vv2,T2) - VT_enthalpy(model,vl2,T2))
        #log(p/p2) = (ΔHvap/R̄)(1/T2 - 1/T)
        #log(p/p2)*R̄/ΔHvap = 1/T - 1/T2
        T3 = 1/(log(p2/p)*R̄/ΔHvap + 1/T2)
        (_,vl3,vv3) = saturation_pressure(model,T2,ChemPotVSaturation(crit_retry = false))
        res,converged = try_sat_temp(model,p,T3,vl3,vv3,scales,method)
        converged && return res
    end
    #no luck here, we need to calculate the critical point
    if isnothing(crit)
        crit = crit_pure(model)
    end
    Tc,pc,vc = crit
    p > pc && return fail
    T2 >= Tc && return fail
    
    if 0.999pc > p > 0.95pc 
        #you could still perform another iteration from a better initial point
        Vl2,Vv2 = x0_sat_pure_crit(model,0.99T2,Tc,pc,vc)
        res,converged = try_sat_temp(model,p,0.99T2,Vl2,Vv2,scales,method)
        converged && return res
    elseif p < 0.5pc
        #very low pressure, we need a better aproximation. luckily we now have a better T
        #TODO: look for a way to detect this case without looking at the critical point
        Vl2,Vv2 = x0_sat_pure(model,T2)
        res,converged = try_sat_temp(model,p,T2,Vl2,Vv2,scales,method)
        converged && return res
    end
    
    if p > 0.999pc
    #almost no hope, only ClapeyronSat works in this range. 
        crit_satmethod = ClapeyronSaturation(;crit)
        return saturation_temperature(model,p,crit_satmethod)
    end
    
    return fail
end

function try_sat_temp(model,p,T0,Vl,Vv,scales,method::AntoineSaturation)
    if T0 isa Base.IEEEFloat # MVector does not work on non bits types, like BigFloat
        v0 = MVector((T0,log(Vl),log(Vv)))
    else
        v0 = SizedVector{3,typeof(T0)}((T0,log(Vl),log(Vv)))
    end
    #we solve volumes in a log scale
    f!(F,x) = Obj_Sat_Temp(model,F,x[1],exp(x[2]),exp(x[3]),p,scales,method)
    r = Solvers.nlsolve(f!,v0, LineSearch(Newton()),NEqOptions(method),ForwardDiff.Chunk{3}())
    sol = Solvers.x_sol(r)
    T = sol[1]
    Vl = exp(sol[2])
    Vv = exp(sol[3])
    converged = check_valid_sat_pure(model,p,Vl,Vv,T)
    return (T,Vl,Vv),converged
end


export AntoineSaturation