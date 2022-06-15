struct AntoineSaturation{T} <: SaturationMethod
    Temp::Union{Nothing,T}
    vl::Union{Nothing,T}
    vv::Union{Nothing,T}
end

AntoineSaturation()=AntoineSaturation{Nothing}(nothing,nothing,nothing)
#if a number is provided as initial point, it will instead proceed to solve directly
function saturation_temperature(model::EoSModel, P, T0::Number)
    sat = x0_sat_pure(model,T0)
    return saturation_temperature_impl(model,p,AntoineSaturation(T0,sat[1],sat[3]))
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

function x0_sat_temp(model::EoSModel,p,method::AntoineSaturation)
    A,B,C = antoine_coef(model)
    lnp̄ = log(p / p_scale(model))
    T0 = T_scale(model)*(B/(A-lnp̄)-C)
    Vl,Vv = x0_sat_pure(model,T0)
    return [T0,Vl,Vv]
end

function saturation_temperature_impl(model,p,method::AntoineSaturation)
    scales = scale_sat_pure(model)
    if isnothing(method.Temp)
        v0 = x0_sat_temp(model,p,method)
    elseif isnothing(method.vl) && isnothing(method.vv)
        vl,vv = x0_sat_pure(model,method.T)
        v0 = [method.Temp,vl,vv]
    else
        v0 = [method.Temp,method.vl,method.vv]
    end

    F = zeros(eltype(v0),length(v0))

    f!(F,x) = Obj_Sat_Temp(model,F,x[1],10^x[2],10^x[3],p,scales,method)
    r = Solvers.nlsolve(f!,v0, LineSearch(Newton()))
    sol = Solvers.x_sol(r)
    T = sol[1]
    Vl = 10^sol[2]
    Vv = 10^sol[3]
    return (T,Vl,Vv)
end

export AntoineSaturation