#for use in models that have activity coefficient defined.
function recombine_impl!(model::ActivityModel)
    recombine!(model.puremodel)
    return model
end

function excess_gibbs_free_energy(model::ActivityModel,p,T,z)
    γ = activity_coefficient(model,p,T,z)
    return sum(z[i]*R̄*T*log(γ[i]) for i ∈ @comps)
end

function test_excess_gibbs_free_energy(model::ActivityModel,p,T,z)
    γ = activity_coefficient(model,p,T,z)
    return sum(z[i]*R̄*T*log(γ[i]) for i ∈ @comps)
end

function volume_impl(model::ActivityModel, p, T, z, phase=:unknown, threaded=false, vol0=nothing)
    if hasfield(typeof(model),:puremodel)
        return volume(model.puremodel.model, p, T, z, phase=phase, threaded=threaded, vol0=vol0)
    else
        return volume(BasicIdeal(), p, T, z, phase=phase, threaded=threaded, vol0=vol0)
    end
end
#for use in models that have gibbs free energy defined.
function activity_coefficient(model::ActivityModel,p,T,z)
    X = gradient_type(p,T,z)
    return exp.(Solvers.gradient(x->excess_gibbs_free_energy(model,p,T,x),z)/(R̄*T))::X
end

function test_activity_coefficient(model::ActivityModel,p,T,z)
    X = gradient_type(p,T,z)
    return exp.(Solvers.gradient(x->excess_gibbs_free_energy(model,p,T,x),z)/(R̄*T))::X
end

x0_sat_pure(model::ActivityModel,T) = x0_sat_pure(__act_to_gammaphi(model,x0_sat_pure),T)

function saturation_pressure(model::ActivityModel,T::Real,method::SaturationMethod)
    return saturation_pressure(__act_to_gammaphi(model,saturation_pressure),T,method)
end

function saturation_temperature(model::ActivityModel,T::Real,method::SaturationMethod)
    return saturation_temperature(__act_to_gammaphi(model,saturation_temperature),T,method)
end

function init_preferred_method(method::typeof(saturation_pressure),model::ActivityModel,kwargs)
    return init_preferred_method(method,__act_to_gammaphi(model,method),kwargs)
end

function init_preferred_method(method::typeof(saturation_temperature),model::ActivityModel,kwargs)
    return init_preferred_method(method,__act_to_gammaphi(model,method),kwargs)
end

function idealmodel(model::T) where T <: ActivityModel
    if hasfield(T,:puremodel)
        puremodel = model.puremodel.model
        return idealmodel(model.puremodel.model)
    else
        return BasicIdeal()
    end
end

function a_res(model::ActivityModel,V,T,z)
    Σz = sum(z)
    pures = model.puremodel
    g_pure_res = sum(z[i]*VT_gibbs_free_energy_res(pures[i],V/Σz,T) for i ∈ @comps)
    p = sum(z[i]*pressure(pures[i],V,T) for i ∈ @comps)/Σz
    g_E = excess_gibbs_free_energy(model,p,T,z)
    p_res = p - Σz*R̄*T/V
    return (g_E+g_pure_res-p_res*V)/(Σz*Rgas(model)*T)
end

function mixing(model::ActivityModel,p,T,z,::typeof(enthalpy))
    f(x) = excess_gibbs_free_energy(model,p,x,z)/x
    dfT = Solvers.derivative(f,T)
    return -dfT*T^2
end

function mixing(model::ActivityModel,p,T,z,::typeof(gibbs_free_energy))
    x = z./sum(z)
    return excess_gibbs_free_energy(model,p,T,z)+dot(z,log.(x))*R̄*T
end

function mixing(model::ActivityModel,p,T,z,::typeof(entropy))
    f(x) = excess_gibbs_free_energy(model,p,x,z)/x
    g,dg = Solvers.f∂f(f,T)
    return -dg*T-g
end

function gibbs_solvation(model::ActivityModel,T)
    binary_component_check(gibbs_solvation,model)
    return gibbs_solvation(__act_to_gammaphi(model,gibbs_solvation),T)
end

function lb_volume(model::ActivityModel,z = SA[1.0])
    b = sum(lb_volume(model.puremodel[i])*z[i] for i in @comps)
    return b
end

function T_scale(model::ActivityModel,z=SA[1.0])
    prod(T_scale(model.puremodel[i])^1/z[i] for i in @comps)^(sum(z))
end

function p_scale(model::ActivityModel,z=SA[1.0])
    0.33*R̄*T_scale(model,z)/lb_volume(model,z)
end

function x0_volume_liquid(model::ActivityModel,T,z)
    pures = model.puremodel
    return sum(z[i]*x0_volume_liquid(pures[i],T,SA[1.0]) for i ∈ @comps)
end

function γdγdn(model::ActivityModel,p,T,z)
    storage = DiffResults.JacobianResult(z)
    γ(_z) = activity_coefficient(model,p,T,_z)
    ForwardDiff.jacobian!(storage,γ,z)
    γz = DiffResults.value(storage)
    dyz = DiffResults.jacobian(storage)
    return γz,dyz
end

#convert ActivityModel into a RestrictedEquilibriaModel
function __act_to_gammaphi(model::ActivityModel,method,ignore = false)
    components = model.components
    if hasfield(typeof(model),:puremodel) && !ignore && model.puremodel.model isa IdealModel
        ActivitySaturationError(model,method)
    end

    if hasfield(typeof(model),:puremodel)
        pure = model.puremodel
        if pure.model isa CompositeModel
            pure = EoSVectorParam(pure.model.fluid,model.components)
        end
    else
        if ignore
            pure = EoSVectorParam(BasicIdeal(),model.components)
        else
            ActivitySaturationError(model,method)
        end
    end
    γϕmodel = GammaPhi(components,model,pure)
end

function bubble_pressure(model::ActivityModel,T,x,method::BubblePointMethod)
    compmodel = __act_to_gammaphi(model,method)
    return bubble_pressure(compmodel,T,x,method)
end

function bubble_temperature(model::ActivityModel,p,x,method::BubblePointMethod)
    compmodel = __act_to_gammaphi(model,method)
    return bubble_temperature(compmodel,p,x,method)
end

function dew_pressure(model::ActivityModel,T,y,method::DewPointMethod)
    compmodel = __act_to_gammaphi(model,method)
    return dew_pressure(compmodel,T,y,method)
end

function dew_temperature(model::ActivityModel,p,y,method::DewPointMethod)
    compmodel = __act_to_gammaphi(model,method)
    return dew_temperature(compmodel,p,y,method)
end

function init_preferred_method(method::typeof(bubble_pressure),model::ActivityModel,kwargs)
    gas_fug = get(kwargs,:gas_fug,false)
    poynting = get(kwargs,:poynting,false)
    return ActivityBubblePressure(;gas_fug,poynting,kwargs...)
end

function init_preferred_method(method::typeof(bubble_temperature),model::ActivityModel,kwargs)
    gas_fug = get(kwargs,:gas_fug,false)
    poynting = get(kwargs,:poynting,false)
    return ActivityBubbleTemperature(;gas_fug,poynting,kwargs...)
end

function init_preferred_method(method::typeof(dew_pressure),model::ActivityModel,kwargs)
    gas_fug = get(kwargs,:gas_fug,false)
    poynting = get(kwargs,:poynting,false)
    return ActivityDewPressure(;gas_fug,poynting,kwargs...)
end

function init_preferred_method(method::typeof(dew_temperature),model::ActivityModel,kwargs)
    gas_fug = get(kwargs,:gas_fug,false)
    poynting = get(kwargs,:poynting,false)
    return ActivityDewTemperature(;gas_fug,poynting,kwargs...)
end

function init_preferred_method(method::typeof(tp_flash),model::ActivityModel,kwargs)
    return RRTPFlash(;kwargs...)
end

function __tpflash_cache_model(model::ActivityModel,p,T,z,equilibrium)
    ignore = is_lle(equilibrium)
    compmodel = __act_to_gammaphi(model,tp_flash,ignore)
    PTFlashWrapper(compmodel,p,T,equilibrium)
end

#LLE point. it does not require an imput concentration, because it assumes that activities are pressure-independent.

function LLE(model::ActivityModel,T;v0=nothing)
    if v0 === nothing
        if length(model) == 2
        v0 = [0.25,0.75]
        else
            throw(error("unable to provide an initial point for LLE pressure"))
        end
    end
    len = length(v0)
    Fcache = zeros(eltype(v0),len)
    f!(F,z) = Obj_LLE(model, F, T, z[1], z[2])
    r  = Solvers.nlsolve(f!,v0,LineSearch(Newton()))
    sol = Solvers.x_sol(r)
    x = sol[1]
    xx = sol[2]
    return x,xx
end

function Obj_LLE(model::ActivityModel, F, T, x, xx)
    x = Fractions.FractionVector(x)
    xx = Fractions.FractionVector(xx)
    γₐ = activity_coefficient(model,1e-3,T,x)
    γᵦ = activity_coefficient(model,1e-3,T,xx)
    F .= γᵦ.*xx .- γₐ.*x
    return F
end

export LLE