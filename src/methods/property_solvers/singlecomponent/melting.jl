function obj_melting_pressure(model::CompositeModel,F,T,vs,vl,p̄,T̄)
    μs = VT_chemical_potential(model.solid, vs, T)[1]
    μl = VT_chemical_potential(model.fluid, vl, T)[1]
    ps = pressure(model.solid, vs, T)
    pl = pressure(model.fluid, vl, T)
    F[1] = (μs - μl)/R̄/T̄
    F[2] = log((ps - pl)/p̄)
    return F
end

struct ChemPotMeltingPressure{V} <: ThermodynamicMethod
    v0::V
    check_triple::Bool
    f_limit::Float64,
    atol::Float64,
    rtol::Float64,
    max_iters::Int
end

function ChemPotMeltingPressure(;v0 = nothing,
                                    check_triple = false,
                                    f_limit = 0.0,
                                    atol = 1e-8,
                                    rtol = 1e-12,
                                    max_iters = 10000)

    return ChemPotMeltingPressure(v0,check_triple,f_limit,atol,rtol,max_iters)
end

function init_preferred_method(method::typeof(melting_pressure),model::CompositeModel{<:EoSModel,<:EoSModel},kwargs)
    ChemPotMeltingPressure(;kwargs...)
end


"""
    pm,vs,vl = melting_pressure(model::CompositeModel,T;v0=x0_melting_pressure(model,T))

Calculates the melting pressure of a `CompositeModel` containing a solid and fluid phase EoS, at a specified pressure.
You can pass a tuple of initial values for the volumes `(vs0,vl0)`.

returns:
- Melting Pressure [`Pa`]
- melting solid volume at specified temperature [`m³`]
- melting liquid volume at specified temperature [`m³`]
"""
function melting_pressure(model::CompositeModel,T,kwargs...)
    method = init_preferred_method(melting_pressure,model,kwargs)
    return melting_pressure(model,T,method)
end

function melting_pressure(model::CompositeModel,T,method::ThermodynamicMethod)
    T = T*T/T
    return melting_pressure_impl(model,T,method)
end
function melting_pressure_impl(model::CompositeModel,T,method::ChemPotMeltingPressure)
    T̄ = T_scale(model.fluid)
    p̄ = p_scale(model.fluid)
    if method.v0 == nothing
        v0 = x0_sublimation_pressure(model,T)
    else
        v0 = method.v0
    end
    V0 = vec2(log(v0[1]),log(v0[2]),T)
    f!(F,x) = obj_melting_pressure(model,F,T,exp(x[1]),exp(x[2]),p̄,T̄)
    results = Solvers.nlsolve(f!,V0)
    x = Solvers.x_sol(results)
    vs = exp(x[1])
    vl = exp(x[2])
    pfus = pressure(model.fluid, vl, T)
    return pfus, vs, vl
end

function x0_melting_pressure(model::CompositeModel,T)
    vs0 = lb_volume(model.solid)*6*1.05/π
    vl0 = lb_volume(model.fluid)*6*1.25/π
    return (vs0,vl0)
end