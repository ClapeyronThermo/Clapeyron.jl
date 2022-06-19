struct AntoineSaturation{T} <: SaturationMethod
    Temp::Union{Nothing,T}
    vl::Union{Nothing,T}
    vv::Union{Nothing,T}
end

AntoineSaturation()=AntoineSaturation{Nothing}(nothing,nothing,nothing)
#if a number is provided as initial point, it will instead proceed to solve directly
function saturation_temperature(model::EoSModel, p, T0::Number)
    sat = x0_sat_pure(model,T0) .|> exp10
    return saturation_temperature_impl(model,p,AntoineSaturation(T0,sat[1],sat[2]))
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
    Vl,Vv = x0_sat_pure(model,T0) .|> exp10
    return (T0,Vl,Vv)
end

#in case that there isn't any antoine coefficients:
#We aproximate to RK, use the cubic antoine, and perform refinement with one Clapeyron Saturation iteration 
function x0_saturation_temperature(model::EoSModel,p,::Nothing)
    Tc,Pc,Vc = crit_pure(model)
    @show Tc
    A,B,C = (6.668322465137264,6.098791871032391,-0.08318016317721941)
    if Pc < p
        nan = zero(p)/zero(p)
        return (nan,nan,nan)
    end
    lnp̄ = log(p / Pc)
    T0 = Tc*(B/(A-lnp̄)-C)
    pii,vli,vvi = saturation_pressure(model,T0)
    
    if isnan(pii)
        nan = zero(p)/zero(p)
        return (nan,nan,nan)
    end

    Δp = (p-pii)
    S_v = VT_entropy(model,vvi,T0)
    S_l = VT_entropy(model,vli,T0)
    ΔS = S_v - S_l
    ΔV = vvi - vli
    dpdt = ΔS/ΔV #≈ (p - pii)/(T-Tnew)
    T = T0 + Δp/dpdt
    vv = volume_virial(model,p,T)
    vl = 0.3*lb_volume(model) + 0.7*vli
    return (T,vl,vv)
end

function saturation_temperature_impl(model,p,method::AntoineSaturation)    
    scales = scale_sat_pure(model)
    if isnothing(method.Temp)
        T0,Vl,Vv = x0_saturation_temperature(model,p)
        Vl,Vv = log10(Vl),log10(Vv) 
    elseif isnothing(method.vl) && isnothing(method.vv)
        Vl,Vv = x0_sat_pure(model,method.T) #exp10
        T0 = method.Temp
    else
        T0,Vl,Vv = method.Temp,method.vl,method.vv
        Vl,Vv = log10(Vl),log10(Vv) 
    end

    T0,Vl,Vv = promote(T0,Vl,Vv)
    if isnan(T0)
        return (T0,T0,T0)
    end
    if T0 isa Base.IEEEFloat # MVector does not work on non bits types, like BigFloat
        v0 = MVector((T0,Vl,Vv))
    else
        v0 = SizedVector{3,typeof(T0)}((T0,Vl,Vv))
    end

    f!(F,x) = Obj_Sat_Temp(model,F,x[1],exp10(x[2]),exp10(x[3]),p,scales,method)
    r = Solvers.nlsolve(f!,v0, LineSearch(Newton()))
    sol = Solvers.x_sol(r)
    T = sol[1]
    Vl = exp10(sol[2])
    Vv = exp10(sol[3])
    valid = check_valid_sat_pure(model,p,Vl,Vv,T)
    if valid
        return (T,Vl,Vv)
    else
        return saturation_temperature_impl(model,p,ClapeyronSaturation())
    end
end

export AntoineSaturation