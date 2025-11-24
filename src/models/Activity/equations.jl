#for use in models that have activity coefficient defined.
function recombine_impl!(model::ActivityModel)
    recombine!(model.puremodel)
    return model
end

function lnγ_impl! end

has_lnγ_impl(model::T) where T = hasmethod(lnγ_impl!,Tuple{Any,T,Any,Any,Any})

function excess_gibbs_free_energy(model::ActivityModel,p,T,z)
    if has_lnγ_impl(model)
        lnγx = lnγ(model,p,T,z)
        return Rgas(model)*T*dot(z,lnγx)
    else
        γ = activity_coefficient(model,p,T,z)
        return Rgas(model)*T*sum(z[i]*log(γ[i]) for i ∈ @comps)
    end
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
    lnγx = lnγ(model,p,T,z)
    if ismutable(lnγx)
        lnγx .= exp.(lnγx)
        return lnγx
    else
        return exp.(lnγ)
    end
end

__γ_unwrap(model) = model

@newmodelsingleton IdealLiquidSolution ActivityModel
excess_gibbs_free_energy(::IdealLiquidSolution,p,T,z) = zero(Base.promote_eltype(T,z))
function lnγ_impl!(res,::IdealLiquidSolution,p,T,z)
    res .= 0
    return res
end

K0_lle_init(::IdealLiquidSolution,p,T,z) = throw(error("IdealLiquidSolution() does not support LLE equilibria."))

function lnγ(model::ActivityModel,p,T,z,cache::TT = nothing) where TT
    X = gradient_type(model,T,z)
    nc = length(z)
    if has_lnγ_impl(model)
        if cache isa Tuple
            result,aux,lnγ,∂lnγ∂n,∂lnγ∂T,_,_,hconfig = cache
            lnγ_impl!(lnγ,model,p,T,z)
            return lnγ
        elseif cache isa AbstractVector
            lnγ_impl!(cache,model,p,T,z)
            return cache
        else
            out = similar(X,nc)
            lnγ_impl!(out,model,p,T,z)
            return out
        end
    else
        fun(x) = excess_gibbs_free_energy(model,p,T,@view(x[1:nc]))/(Rgas(model)*T)
        if cache isa Tuple
            result,aux,lnγ,∂lnγ∂n,∂lnγ∂T,_,_,hconfig = cache
            aux .= 0
            aux[1:nc] = z
            gconfig = Solvers._GradientConfig(hconfig)
            _result = ForwardDiff.gradient!(result, fun, aux, gconfig, Val{false}())
            dresult = DiffResults.gradient(_result)
            lnγx = @view dresult[1:nc]
            lnγ .= lnγx
            return lnγ
        elseif cache isa AbstractVector
            ForwardDiff.gradient!(cache,fun,z)
            return cache
        else
            return Solvers.gradient(fun,z)
        end
    end
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

function ∂lnγ∂n(model,p,T,z,cache = nothing)
    nc = length(z)
    RT = Rgas(model)*T
    n = sum(z)
    fun_g(w) = excess_gibbs_free_energy(model,p,T,@view(w[1:nc]))/RT
    function fun_lnγ(out,w)
        Clapeyron.lnγ(model,p,T,@view(w[1:nc]),@view(out[1:nc]))
        return out
    end
    if cache == nothing
        if has_lnγ_impl(model)
            lnγ = zeros(Base.promote_eltype(model,p,T,z),nc)
            ∂lnγ∂ni = ForwardDiff.jacobian!(lnγ,fun_lnγ,z)
            g_E = dot(z,lnγ)*RT
            return g_E,lnγ,∂lnγ∂ni
        else
            hresult = DiffResults.HessianResult(z)
            result = ForwardDiff.hessian!(hresult,fun_g,z)
            g_E = DiffResults.value(result)*RT
            lnγ = DiffResults.gradient(result)
            ∂lnγ∂ni = DiffResults.hessian(result)
            return g_E,lnγ,∂lnγ∂ni
        end
    else
        result,aux,lnγ,∂lnγ∂ni,_,_,_,hconfig,jcache = cache
        aux .= 0
        aux[1:nc] .= z
        if has_lnγ_impl(model)
            jconfig = Solvers._JacobianConfig(hconfig)
            jresult = ForwardDiff.DiffResults.MutableDiffResult(result.derivs[1],(result.derivs[2],))
            _result = ForwardDiff.jacobian!(jresult,fun_lnγ,jcache,aux,jconfig,Val{false}())
            ∂lnγ = DiffResults.jacobian(_result)
            ∂lnγ∂ni .=  @view ∂lnγ[1:nc,1:nc]
            ∂g_E = DiffResults.value(_result)
            lnγ .= @view ∂g_E[1:nc]
            g_E = dot(z,lnγ)*RT
            return g_E,lnγ,∂lnγ∂ni
        else
            _result = ForwardDiff.hessian!(result, fun_g, aux, hconfig, Val{false}())
            g_E = DiffResults.value(_result)*RT
            ∂g_E = DiffResults.gradient(_result)
            lnγ .= @view ∂g_E[1:nc]
            ∂lnγ = DiffResults.hessian(_result)
            ∂lnγ∂ni .= @view ∂lnγ[1:nc,1:nc]
            return g_E,lnγ,∂lnγ∂ni
        end
    end
end

function ∂lnγ∂n∂T(model,p,T,z,cache = nothing)
    nc = length(z)
    RT = Rgas(model)*T
    n = sum(z)
    fun_g(w) = excess_gibbs_free_energy(model,p,w[nc+1],@view(w[1:nc]))/RT
    function fun_lnγ(out,w)
        Clapeyron.lnγ(model,p,w[1:nc+1],@view(w[1:nc]),@view(out[1:nc]))
        return out
    end
    if cache == nothing
        if has_lnγ_impl(model)
            lnγ = zeros(Base.promote_eltype(model,p,T,z),nc)
            aux = similar(lnγ,nc+1)
            aux[1:nc] = z
            aux[nc+1] = T
            ∂g_E = ForwardDiff.jacobian!(lnγ,fun_lnγ,aux)
            ∂lnγ∂ni = ∂g_E[1:nc,1:nc]
            ∂lnγ∂T = resize!(aux,nc)
            ∂lnγ∂T .= @view ∂g_E[:,nc + 1]
            g_E = dot(z,lnγ)*RT
            return g_E,lnγ,∂lnγ∂ni,∂lnγ∂T
        else
            hresult = DiffResults.HessianResult(z)
            result = ForwardDiff.hessian!(hresult,fun_g,z)

            g_E = DiffResults.value(result)*RT
            lnγ_and_T = DiffResults.gradient(result)
            lnγ = resize!(lnγ_and_T,nc)

            ∂lnγ∂ni∂T = DiffResults.hessian(result)
            ∂lnγ∂ni = ∂lnγ∂ni∂T[1:nc,1:nc]
            ∂lnγ∂T = ∂lnγ∂ni∂T[:,nc + 1]
            return g_E,lnγ,∂lnγ∂ni,∂lnγ∂T
        end
    else
        result,aux,lnγ,∂lnγ∂ni,∂lnγ∂T,_,_,hconfig,jcache = cache
        aux .= 0
        aux[1:nc] .= z
        aux[nc+1] = T
        if has_lnγ_impl(model)
            jconfig = Solvers._JacobianConfig(hconfig)
            jresult = ForwardDiff.DiffResults.MutableDiffResult(result.derivs[1],(result.derivs[2],))
            _result = ForwardDiff.jacobian!(jresult,fun_lnγ,jcache,aux,jconfig,Val{false}())
            ∂lnγ = DiffResults.jacobian(_result)
            ∂lnγ∂ni .=  @view ∂lnγ[1:nc,1:nc]
            ∂lnγ∂T .= @view ∂lnγ[:,nc + 1]
            ∂g_E = DiffResults.value(_result)
            lnγ .= @view ∂g_E[1:nc]
            g_E = dot(z,lnγ)*RT
            return g_E,lnγ,∂lnγ∂ni,∂lnγ∂T
        else
            _result = ForwardDiff.hessian!(result, fun_g, aux, hconfig, Val{false}())
            g_E = DiffResults.value(_result)*RT
            ∂g_E = DiffResults.gradient(_result)
            lnγ .= @view ∂g_E[1:nc]
            ∂lnγ = DiffResults.hessian(_result)
            ∂lnγ∂ni .= @view ∂lnγ[1:nc,1:nc]
            ∂lnγ∂T .= @view ∂lnγ[1:nc,nc + 1]
            return g_E,lnγ,∂lnγ∂ni,∂lnγ∂T
        end
    end
end


__act_to_gammaphi(model::ActivityModel) = __act_to_gammaphi(model,nothing,true)
GammaPhi(model::ActivityModel) = __act_to_gammaphi(model)
#convert ActivityModel into a RestrictedEquilibriaModel
function __act_to_gammaphi(model::ActivityModel,method,ignore = false)
    components = component_list(model)
    if hasfield(typeof(model),:puremodel) && !ignore && model.puremodel.model isa IdealModel
        ActivitySaturationError(model,method)
    end

    if hasfield(typeof(model),:puremodel)
        pure = model.puremodel
        if pure.model isa CompositeModel
            pure = EoSVectorParam(pure.model.fluid,components)
        end
    else
        if ignore
            pure = EoSVectorParam(BasicIdeal(),components)
        else
            ActivitySaturationError(model,method)
        end
    end
    γϕmodel = GammaPhi(components,model,pure)
end

for f in (:bubble_pressure,:bubble_temperature,:dew_pressure,:dew_temperature)
    @eval begin
        function $f(model::ActivityModel,T,x,method::ThermodynamicMethod)
            compmodel = __act_to_gammaphi(model,method)
            return $f(compmodel,T,x,method)
        end
    end
end

function init_preferred_method(method::typeof(bubble_pressure),model::ActivityModel,kwargs)
    return ActivityBubblePressure(;kwargs...)
end

function init_preferred_method(method::typeof(bubble_temperature),model::ActivityModel,kwargs)
    return FugBubbleTemperature(;kwargs...)
end

function init_preferred_method(method::typeof(dew_pressure),model::ActivityModel,kwargs)
    return ActivityDewPressure(;kwargs...)
end

function init_preferred_method(method::typeof(dew_temperature),model::ActivityModel,kwargs)
    return FugDewTemperature(;kwargs...)
end

function init_preferred_method(method::typeof(tp_flash),model::ActivityModel,kwargs)
    return MichelsenTPFlash(;kwargs...)
end

function gas_model(model::T) where T <:ActivityModel
    return gas_model(__act_to_gammaphi(model,tp_flash,true))
end

function PTFlashWrapper(model::ActivityModel,p,T,z,equilibrium)
    ignore = is_lle(equilibrium)
    compmodel = __act_to_gammaphi(model,tp_flash,ignore)
    return PTFlashWrapper(compmodel,p,T,z,equilibrium)
end

function __tpflash_cache_model(model::ActivityModel,p,T,z,equilibrium)
    PTFlashWrapper(model,p,T,z,equilibrium)
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

function PT_property(model::ActivityModel,p,T,z,phase,threaded,vol0,f::F,v::Val{UseP}) where {F,UseP}
    γϕ = __act_to_gammaphi(model)
    PT_property(γϕ,p,T,z,phase,threaded,vol0,f,v)
end

function set_reference_state!(model::ActivityModel,reference_state::ReferenceState;verbose = verbose)
    γϕ = __act_to_gammaphi(model)
    set_reference_state!(γϕ,reference_state;verbose)
end

reference_state(model::ActivityModel) = reference_state(model.puremodel)