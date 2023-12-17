function obj_sublimation_pressure(model::CompositeModel,F,T,vs,vv,p̄,T̄)
    z = SA[1.0]
    eos_solid(V) = eos(model.solid,V,T,z)
    eos_fluid(V) = eos(model.fluid,V,T,z)
    A_v,Av_v = Solvers.f∂f(eos_fluid,vv)
    A_s,Av_s =Solvers.f∂f(eos_solid,vs)
    μv = muladd(-vv,Av_v,A_v)
    μs = muladd(-vs,Av_s,A_s)
    ps = - Av_s
    pv = - Av_v
    #=
    μs = VT_chemical_potential(model.solid, vs, T)[1]
    μv = VT_chemical_potential(model.fluid, vv, T)[1]
    ps = pressure(model.solid, vs, T)
    pv = pressure(model.fluid, vv, T) =#
    F[1] = (μs - μv)/R̄/T̄
    F[2] = (ps - pv)/p̄
    return F
end

struct ChemPotSublimationPressure{V} <: ThermodynamicMethod
    v0::V
    check_triple::Bool
    f_limit::Float64,
    atol::Float64,
    rtol::Float64,
    max_iters::Int
end

function ChemPotSublimationPressure(;v0 = nothing,
                                    check_triple = false,
                                    f_limit = 0.0,
                                    atol = 1e-8,
                                    rtol = 1e-12,
                                    max_iters = 10000)

    return ChemPotSublimationPressure(v0,check_triple,f_limit,atol,rtol,max_iters)
end

function init_preferred_method(method::typeof(sublimation_pressure),model::CompositeModel{<:EoSModel,<:EoSModel},kwargs)
    ChemPotSublimationPressure(;kwargs...)
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
    return sublimation_pressure(model,p,method)
end

function sublimation_pressure(model,T,method::ThermodynamicMethod)
    single_component_check(sublimation_pressure,model)
    T = T*T/T
    return sublimation_pressure_impl(model,T,method)
end

function sublimation_pressure_impl(model::CompositeModel,T,method::ChemPotSublimationPressure)
    T̄ = T_scale(model.fluid)
    p̄ = p_scale(model.fluid)
    if method.v0 == nothing
        v0 = x0_sublimation_pressure(model,T)
    else
        v0 = method.v0
    end
    V0 = vec2(log(v0[1]),log(v0[2]),T)
    f!(F,x) = obj_sublimation_pressure(model,F,T,exp10(x[1]),exp10(x[2]),p̄,T̄)
    results = Solvers.nlsolve(f!,V0)
    x = Solvers.x_sol(results)
    vs = exp(x[1])
    vv = exp(x[2])
    psub = pressure(model.fluid, vv, T)
    return psub, vs, vv
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
    vs0 = volume(solid,P0,T,z,vol0 = vs_at_0)
    vv0 = R̄*T/P0
    return (vs0,vv0)
end
