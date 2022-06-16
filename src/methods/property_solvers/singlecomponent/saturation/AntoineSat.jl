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

x0_saturation_temperature(model,p) = x0_sat_temperature(model,p,AntoineSaturation())

function x0_saturation_temperature(model::EoSModel,p,::AntoineSaturation)
    A,B,C = antoine_coef(model)
    lnp̄ = log(p / p_scale(model))
    T0 = T_scale(model)*(B/(A-lnp̄)-C)
    Vl,Vv = x0_sat_pure(model,T0) .|> exp10
    return (T0,Vl,Vv)
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
    return (T,Vl,Vv)
end

export AntoineSaturation