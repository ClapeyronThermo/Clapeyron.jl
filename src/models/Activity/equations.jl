#for use in models that have activity coefficient defined.
function recombine_impl!(model::ActivityModel)
    recombine!(model.puremodel)
    return model
end

function excess_gibbs_free_energy(model::ActivityModel,p,T,z)
    γ = activity_coefficient(model,p,T,z)
    return Rgas(model)*T*sum(z[i]*log(γ[i]) for i ∈ @comps)
end

function test_excess_gibbs_free_energy(model::ActivityModel,p,T,z)
    γ = activity_coefficient(model,p,T,z)
    return Rgas(model)*T*sum(z[i]*log(γ[i]) for i ∈ @comps)
end

function volume_impl(model::ActivityModel, p, T, z, phase, threaded, vol0)
    if hasfield(typeof(model),:puremodel)
        return volume(model.puremodel.model, p, T, z, phase=phase, threaded=threaded, vol0=vol0)
    else
        return volume(BasicIdeal(), p, T, z, phase=phase, threaded=threaded, vol0=vol0)
    end
end
#for use in models that have Gibbs energy defined.
function activity_coefficient(model::ActivityModel,p,T,z)
    X = gradient_type(model,T+p,z)
    return exp.(Solvers.gradient(x->excess_gibbs_free_energy(model,p,T,x),z)/(Rgas(model)*T))::X
end

function activity_coefficient_impl(model::ActivityModel,p,T,z,μ_ref,reference,phase,threaded,vol0)
    #TODO: what to do if the reference is not pure?
    return activity_coefficient(model,p,T,z)
end

reference_chemical_potential_type(model::ActivityModel) = :zero

function activity(model::ActivityModel,p,T,z)
    γ = activity_coefficient(model,p,T,z)
    ∑z = sum(z)
    return γ .* z ./ ∑z
end

function activity_impl(model::ActivityModel,p,T,z,μ_ref,reference,phase,threaded,vol0)
    #TODO: what to do if the reference is not pure?
    return activity(model,p,T,z)
end

function test_activity_coefficient(model::ActivityModel,p,T,z)
    X = gradient_type(model,T+p,z)
    return exp.(Solvers.gradient(x->excess_gibbs_free_energy(model,p,T,x),z)/(R̄*T))::X
end

saturation_model(model::ActivityModel) = __act_to_gammaphi(model,saturation_model)

function idealmodel(model::T) where T <: ActivityModel
    if hasfield(T,:puremodel)
        puremodel = model.puremodel.model
        return idealmodel(model.puremodel.model)
    else
        return BasicIdeal()
    end
end

#=

this is technically wrong on the strict sense of helmholtz residual energy,
but allows us to evaluate the excess terms of an activity model with ease.

The main problem is that activity models are defined in a P-T basis, while the Helmholtz energy framework used by Clapeyron requires a V-T basis.
we circunvent this by using the dispatches on PT_property.
Activity models are transformed into a GammaPhi wrapper that evaluates the pure and excess parts in a correct way.

=#
function eos_impl(model::ActivityModel,V,T,z)
    return excess_gibbs_free_energy(model,V,T,z) + reference_state_eval(model,V,T,z)
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

function lb_volume(model::ActivityModel,z)
    b = sum(lb_volume(model.puremodel[i])*z[i] for i in @comps)
    return b
end

function lb_volume(model::ActivityModel,T,z)
    b = sum(lb_volume(model.puremodel[i],T,SA[1.0])*z[i] for i in @comps)
    return b
end

function T_scale(model::ActivityModel,z)
    prod(T_scale(model.puremodel[i])^1/z[i] for i in @comps)^(sum(z))
end

function p_scale(model::ActivityModel,z)
    T = T_scale(model,z)
    0.33*R̄*T/lb_volume(model,T,z)
end

function x0_volume_liquid(model::ActivityModel,p,T,z)
    pures = model.puremodel
    return sum(z[i]*x0_volume_liquid(pures[i],p,T,SA[1.0]) for i ∈ @comps)
end

function γdγdn(model::ActivityModel,p,T,z)
    storage = DiffResults.JacobianResult(z)
    γ(_z) = activity_coefficient(model,p,T,_z)
    ForwardDiff.jacobian!(storage,γ,z)
    γz = DiffResults.value(storage)
    dyz = DiffResults.jacobian(storage)
    return γz,dyz
end

__act_to_gammaphi(model::ActivityModel) = __act_to_gammaphi(model,nothing,true)
GammaPhi(model::ActivityModel) = __act_to_gammaphi(model)
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

#LLE point. It does not require an input concentration, because it assumes that activities are pressure-independent.
"""
    LLE(model::ActivityModel, T; v0=nothing)

Calculates the Liquid-Liquid equilibrium compositions at a given temperature `T` in `[K]`.

Returns a tuple, containing:
- Liquid composition `x₁`
- Liquid composition `x₂`

`v0` is a vector containing `vcat(x1[1:nc-1],x2[1:nc-1])`.
"""
function LLE(model::ActivityModel,T;v0=nothing)
    nc = length(model)
    vv0 = zeros(Base.promote_eltype(model,T),2*nc-2)
    if v0 === nothing
        if nc == 2
            vv0 .= [0.25,0.75]
        else
            throw(error("unable to provide an initial point for LLE pressure"))
        end
    else
        vv0 = zeros(eltype())
        if 2*length(model) == length(v0)
            vv0[1:nc-1] .= v0[1:nc-1]
            vv0[nc:end] .= v0[(nc+1):(2*nc-1)]
        else
            vv0 .= v0
        end
    end

    len = length(vv0)
    f!(F,z) = Obj_LLE(model, F, T, @view(z[1:nc-1]), @view(z[nc:end]))
    r  = Solvers.nlsolve(f!,vv0,LineSearch(Newton()))
    sol = Solvers.x_sol(r)
    x = FractionVector(sol[1:nc-1]) |> collect
    xx = FractionVector(sol[nc:end]) |> collect
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

function tpd(model::ActivityModel,p,T,z,cache = tpd_cache(model,p,T,z);reduced = false,break_first = false,lle = false,tol_trivial = 1e-5,strategy = :pure, di = nothing)
    #TODO: support tpd with vle and activities?
    if !lle
        throw(ArgumentError("tpd only supports lle search with Activity Models. try using `tpd(model,p,T,z,lle = true)`"))
    end
    γϕmodel = __act_to_gammaphi(model,tpd,true)
    return tpd(γϕmodel,p,T,z,cache;reduced,break_first,lle,tol_trivial,strategy,di)
end

function PT_property(model::ActivityModel,p,T,z,phase,threaded,vol0,f::F,v::Val{UseP}) where {F,UseP}
    γϕ = __act_to_gammaphi(model)
    PT_property(γϕ,p,T,z,phase,threaded,vol0,f,v)
end

function set_reference_state!(model::ActivityModel,reference_state::ReferenceState;verbose = verbose)
    γϕ = __act_to_gammaphi(model)
    set_reference_state!(γϕ,reference_state;verbose)
end

reference_state(model::ActivityModel) = reference_state(model.puremodel)