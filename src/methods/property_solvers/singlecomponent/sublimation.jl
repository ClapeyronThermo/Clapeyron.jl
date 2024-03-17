function obj_sublimation_pressure(model::CompositeModel,T,vs,vv,p_scale,μ_scale)
    solid = solid_model(model)
    fluid = fluid_model(model)
    return μp_equality1_p(solid,fluid,vs,vl,T,p_scale,μ_scale)
end

struct ChemPotSublimationPressure{V} <: ThermodynamicMethod
    v0::V
    check_triple::Bool
    f_limit::Float64
    atol::Float64
    rtol::Float64
    max_iters::Int
end

function ChemPotSublimationPressure(;v0 = nothing,
                                    check_triple = false,
                                    f_limit = 0.0,
                                    atol = 1e-8,
                                    rtol = 1e-12,
                                    max_iters = 100)

    return ChemPotSublimationPressure(v0,check_triple,f_limit,atol,rtol,max_iters)
end

"""
    psub,vs,vv = sublimation_pressure(model::CompositeModel,T;v0=x0_sublimation_pressure(model,T))

Calculates the sublimation pressure of a `CompositeModel` containing a solid and fluid phase EoS, at a specified pressure.
You can pass a tuple of initial values for the volumes `(vs0,vv0)`.

returns:
- Sublimation Pressure [`Pa`]
- Sublimation solid volume at specified temperature [`m³`]
- Sublimation vapour volume at specified temperature [`m³`]
"""
function sublimation_pressure(model::CompositeModel,T;kwargs...)
    method = init_preferred_method(sublimation_pressure,model,kwargs)
    return sublimation_pressure(model,T,method)
end

function init_preferred_method(method::typeof(sublimation_pressure),model::CompositeModel{<:EoSModel,<:EoSModel},kwargs)
    ChemPotSublimationPressure(;kwargs...)
end

function sublimation_pressure(model,T,method::ThermodynamicMethod)
    single_component_check(sublimation_pressure,model)
    T = T*T/T
    return sublimation_pressure_impl(model,T,method)
end

function sublimation_pressure_impl(model::CompositeModel,T,method::ChemPotSublimationPressure)
    if method.v0 == nothing
        v0 = x0_sublimation_pressure(model,T)
    else
        v0 = method.v0
    end
    vs0,vv0 = v0
    _0 = zero(vs0*vv0*T*oneunit(eltype(model)))
    nan = _0/_0
    fail = (nan,nan,nan)
    valid_input = check_valid_2ph_input(vs0,vv0,true,T)
    if !valid_input
        return fail
    end
    fluid = fluid_model(model)
    solid = solid_model(model)
    ps,μs = scale_sat_pure(fluid)
    result,converged = try_2ph_pure_pressure(solid,fluid,T,vs0,vv0,ps,μs,method)
    if converged
        return result
    else
        return fail
    end
end

function x0_sublimation_pressure(model,T)
    #we can suppose we are in a low pressure regime, we treat the solid as a liquid,
    #and apply the zero pressure aproximation.
    solid = solid_model(model)
    fluid = fluid_model(model)
    R̄ = Rgas(solid)
    z = SA[1.0]
    vs_at_0 = volume(solid,0.0,T,phase = :s)
    ares = a_res(solid, vs_at_0, T, z)
    lnϕ_s0 = ares - 1 + log(R̄*T/vs_at_0)
    P0 = exp(lnϕ_s0)
    vv0 = R̄*T/P0
    vs0 = vs_at_0
    return vs0,vv0
end


function Obj_Sub_Temp(model::EoSModel, F, T, V_s, V_v,p,p̄,T̄)
    z = SA[1.0]
    eos_solid(V) = eos(model.solid,V,T,z)
    eos_fluid(V) = eos(model.fluid,V,T,z)
    A_v,Av_v = Solvers.f∂f(eos_fluid,V_v)
    A_s,Av_s =Solvers.f∂f(eos_solid,V_s)
    g_v = muladd(-V_v,Av_v,A_v)
    g_s = muladd(-V_s,Av_s,A_s)

    F1 = -(Av_v+p)/p̄
    F2 = -(Av_s+p)/p̄
    F3 = (g_v-g_s)/(R̄*T̄)
    return SVector(F1,F2,F3)
end

struct ChemPotSublimationTemperature{V} <: ThermodynamicMethod
    T0::Union{Nothing,V}
    v0::V
    check_triple::Bool
    f_limit::Float64
    atol::Float64
    rtol::Float64
    max_iters::Int
end

function ChemPotSublimationTemperature(;v0 = nothing,
                                    T0 = nothing,
                                    check_triple = false,
                                    f_limit = 0.0,
                                    atol = 1e-8,
                                    rtol = 1e-12,
                                    max_iters = 10000)

    return ChemPotSublimationTemperature(v0,T0,check_triple,f_limit,atol,rtol,max_iters)
end

"""
    pm,vs,vl = sublimation_temperature(model::CompositeModel,T;v0=x0_sublimation_pressure(model,T))

Calculates the sublimation temperature of a `CompositeModel` containing a solid and fluid phase EoS, at a specified pressure.
You can pass a tuple of initial values for the volumes `(vs0,vl0)`.

returns:
- Sublimation Temperature [`K`]
- sublimation solid volume at specified pressure [`m³`]
- sublimation vapour volume at specified pressure [`m³`]
"""
function sublimation_temperature(model::CompositeModel,p;kwargs...)
    method = init_preferred_method(sublimation_temperature,model,kwargs)
    return sublimation_temperature(model,p,method)
end
function init_preferred_method(method::typeof(sublimation_temperature),model::CompositeModel{<:EoSModel,<:EoSModel},kwargs)
    ChemPotSublimationTemperature(;kwargs...)
end

function sublimation_temperature(model::CompositeModel,p,method::ThermodynamicMethod)
    p = p*p/p
    return sublimation_temperature_impl(model,p,method)
end

function sublimation_temperature_impl(model::CompositeModel,p,method::ChemPotSublimationTemperature)
    solid = solid_model(model)
    fluid = fluid_model(model)
    T̄ = T_scale(fluid)
    p̄ = p_scale(fluid)
    if method.v0 == nothing
        v0 = x0_sublimation_temperature(model,p)
    else
        v0 = method.v0
    end
    _1 = oneunit(p*1.0*one(eltype(solid))*one(eltype(fluid)))
    V0 = SVector(v0[1]*_1,log(v0[2])*_1,log(v0[3])*_1)
    f!(F,x) = Obj_Sub_Temp(model,F,x[1],exp(x[2]),exp(x[3]),p,p̄,T̄)
    results = Solvers.nlsolve(f!,V0,TrustRegion(Newton(),Dogleg()),NEqOptions(method))
    x = Solvers.x_sol(results)
    vs = exp(x[2])
    vv = exp(x[3])
    Tfus = x[1]
    converged = check_valid_eq2(solid,fluid,p,vs,vv,Tfus)
    if converged
    return Tfus, vs, vv
    else
        nan = zero(Tfus)/zero(Tfus)
        return nan,nan,nan
    end
end

function x0_sublimation_temperature(model::CompositeModel,p)
    trp = triple_point(model)
    pt = trp[2]
    vs0 = trp[3]
    vv0 = trp[5]
    Δv = vv0-vs0
    Tt = trp[1]
    hs0 = VT_enthalpy(model.solid,vs0,Tt)
    hv0 = VT_enthalpy(model.fluid,vv0,Tt)
    Δh = hv0-hs0
    T0 = Tt*exp(Δv*(p-pt)/Δh)
    return T0,vs0,vv0
end
