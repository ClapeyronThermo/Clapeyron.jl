"""
    GammaPhi{γ,Φ} <: RestrictedEquilibriaModel

Wrapper struct to signal that a `CompositeModel` uses an activity model in conjunction with a fluid.
"""
struct GammaPhi{γ,Φ} <: RestrictedEquilibriaModel
    components::Vector{String}
    activity::γ
    fluid::EoSVectorParam{Φ}
end

reference_state(model::GammaPhi) = reference_state(model.fluid)

function Base.show(io::IO,mime::MIME"text/plain",model::GammaPhi)
    print(io,"γ-ϕ Model")
    length(model) == 1 && print(io, " with 1 component:")
    length(model) > 1 && print(io, " with ", length(model), " components:")
    println(io)
    show_pairs(io,model.components)
    act = model.activity
    if hasfield(typeof(act),:puremodel)
        print(io,'\n',"Activity Model: ", parameterless_type(act))
    else
        print(io,'\n',"Activity Model: ",typeof(act))
    end
    print(io,'\n',"Fluid Model: ",typeof(model.fluid.model))
    show_reference_state(io,model;space = true)
end

fluid_model(model::GammaPhi) = model.fluid.model
__γ_unwrap(model::GammaPhi) = __γ_unwrap(model.activity)
@inline gas_model(model::GammaPhi) = gas_model(model.fluid.model)
Base.eltype(model::GammaPhi) = Base.promote_eltype(__γ_unwrap(model),gas_model(model))

function excess_gibbs_free_energy(model::GammaPhi,p,T,z)
    return excess_gibbs_free_energy(model.activity,p,T,z)
end

reference_chemical_potential_type(model::GammaPhi) = reference_chemical_potential_type(model.activity)

function volume_impl(model::GammaPhi,p,T,z,phase,threaded,vol0)
    return volume_impl(model.fluid.model,p,T,z,phase,threaded,vol0)
end

molecular_weight(model::GammaPhi,z) = molecular_weight(model.fluid.model,z)
saturation_model(model::GammaPhi) = saturation_model(model.fluid)
idealmodel(model::GammaPhi) = idealmodel(model.fluid.model)

function init_preferred_method(method::typeof(tp_flash),model::GammaPhi,kwargs)
    MichelsenTPFlash(;kwargs...)
end

# Error handling for Activity models that don't provide saturation properties, in the context of VLE.
function ActivitySaturationError(model,method)
    throw(ArgumentError("$method requires $model to be used along with another EoS model that supports saturation properties. If you are using an Activity Model as a raw input, use `CompositeModel(components, liquid = activity_model, fluid = fluid_model)` instead."))
end

function gibbs_solvation(model::GammaPhi,T)
    z = [1.0,1e-30]
    p,v_l,v_v = saturation_pressure(model.fluid[1],T)
    p2,v_l2,v_v2 = saturation_pressure(model.fluid[2],T)
    γ = activity_coefficient(model,p,T,z)
    K = v_v/v_l*γ[2]*p2/p
    return -R̄*T*log(K)
end

function __calculate_reference_state_consts(model::GammaPhi,v,T,p,z,H0,S0,phase)
    ∑z = sum(z)
    S00 = entropy(model,p,T,z,phase = phase)
    a1 = (S00 - S0)#/∑z
    H00 = enthalpy(model,p,T,z,phase = phase)
    a0 = (-H00 + H0)#/∑z
    return a0,a1
end

function PTFlashWrapper(model::GammaPhi,p,T,z,equilibrium)
    fluidmodel = model.fluid
    #check that we can actually solve the equilibria
    if fluidmodel isa IdealModel && !is_lle(equilibrium)
        ActivitySaturationError(model.activity,tp_flash)
    end
    TT = Base.promote_eltype(model,p,T,z)
    wrapper = PTFlashWrapper{TT}(model,equilibrium,fluidmodel.pure)
    if is_vle(equilibrium) || is_unknown(equilibrium)
        update_temperature!(wrapper,T)
    end
    return wrapper
end

function __tpflash_cache_model(model::GammaPhi,p,T,z,equilibrium)
    PTFlashWrapper(model,p,T,z,equilibrium)
end

function __tpflash_gibbs_reduced(wrapper::PTFlashWrapper{<:GammaPhi},p,T,x,y,β,eq,vols)
    model = wrapper.model
    gibbs = zero(Base.promote_eltype(model,p,T,x,β))
    if !isone(β)
        gx,_ = modified_gibbs(wrapper,p,T,x,:l)
        gibbs += gx*(1-β)
    end

    if is_vle(eq) && !iszero(β)
        vv = vols[2]
        gy,_ = modified_gibbs(wrapper,p,T,y,:v,vols[2])
        gibbs += gy*β
    elseif !iszero(β) #lle
        gy,_ = modified_gibbs(wrapper,p,T,y,:l)
        gibbs += gy*β
    end
    return gibbs
end

function tpd_input_composition(wrapper::PTFlashWrapper{<:GammaPhi},p,T,z,lle,cache = tpd_cache(wrapper,p,T,z,di))

    TT = Base.promote_eltype(wrapper.model,p,T,z)

    pures = wrapper.model.fluid.pure
    model = wrapper.model
    fluidmodel = model.fluid.model
    RT = R̄*T

    d_l,d_v,_,_,_,Hϕ = cache

    TT = Base.promote_eltype(model,p,T,z)

    n = sum(z)
    logsumz = log(n)
    d,vl = tpd_lnϕ_and_v!(last(cache),wrapper,p,T,z,nothing,false,:liquid)
    d_l .= d
    d_l .+= log.(z) .- logsumz

    lle && return copy(d_l),:liquid,vl

    d,vv = tpd_lnϕ_and_v!(last(cache),wrapper,p,T,z,nothing,false,:vapour)
    d_v .= d
    d_v .+= log.(z) .- logsumz
    gr_l = dot(z,d_l)
    gr_v = dot(z,d_v)
    if gr_l < gr_v
        return copy(d_l),:liquid,vl
    else
        return copy(d_v),:vapour,vv
    end
end

function tpd_lnϕ_and_v!(cache,wrapper::PTFlashWrapper,p,T,w,vol0,liquid_overpressure = false,phase = :l,_vol = nothing)
    model = wrapper.model
    RT = R̄*T
    if is_liquid(phase)
        γmodel = __γ_unwrap(model)
        #=
        If the model is not an activity model, then PTFlashWrapper is wrapping
        a normal helmholtz model, we just return lnϕ.
        =#
        if γmodel isa ActivityModel
            logγx = lnγ(γmodel,p,T,w,cache)
            v = zero(eltype(logγx))
            return logγx,v,true
        elseif is_vle(wrapper.equilibrium)
            logγx,v = __lnγ_sat(wrapper,p,T,w,cache)
            return logγx,v,true
        end
    end
    fxy,v,overpressure = tpd_lnϕ_and_v!(cache,gas_model(model),p,T,w,vol0,liquid_overpressure,phase,_vol)
    is_vapour(phase) && !is_lle(wrapper.equilibrium) && tpd_delta_d_vapour!(fxy,wrapper,p,T)

    return fxy,v,overpressure
end

function __lnγ_sat(wrapper::PTFlashWrapper,p,T,w,cache = nothing,vol0 = nothing,vol = volume(wrapper.model,p,T,w,vol0 = vol0,phase = :l))
    model = wrapper.model
    μmix_temp = VT_chemical_potential_res!(cache,model,vol,T,w)
    result,aux,logγ,A1,μmix,x2,x3,hconfig = cache
    μmix .= μmix_temp
    sat = wrapper.sat
    fug = wrapper.fug
    RT = Rgas(model)*T
    for i in 1:length(logγ)
        ϕᵢ = fug[i]
        pᵢ,vpureᵢ,_ = sat[i]

        μᵢ_over_RT = log(ϕᵢ) + log(pᵢ*vpureᵢ/RT)
        logγ[i] = log(vpureᵢ/vol) + μmix[i]/RT - μᵢ_over_RT -  vpureᵢ*(p - pᵢ)/RT
    end
    return logγ,vol
end

function modified_∂lnϕ∂n(wrapper::PTFlashWrapper{<:GammaPhi}, p, T, z, cache; phase = :unknown, vol0 = nothing)
    model = wrapper.model
    if is_vapour(phase)
        lnϕ,∂lnϕ∂n,vol =  modified_∂lnϕ∂n(gas_model(model),p,T,z,cache;phase,vol0)
        tpd_delta_d_vapour!(lnϕ,wrapper,p,T)
        return lnϕ,∂lnϕ∂n,vol
    elseif is_liquid(phase)
        g_E,lnγ,∂lnγ∂ni = ∂lnγ∂n(__γ_unwrap(model),p,T,z,cache)
        return lnγ,∂lnγ∂ni,zero(g_E)
    else
        throw(error("invalid specification for phase: $phase"))
    end
end

function PT_property(model::GammaPhi,p,T,z,phase,threaded,vol0,f::F,USEP::Val{UseP}) where {F,UseP}
    if is_vapour(phase)
        return PT_property(gas_model(model),p,T,z,phase,threaded,vol0,f,USEP)
    else #liquid or unknown
        wrapper = PTFlashWrapper(model,p,T,z,:vle)
        return PT_property(wrapper,p,T,z,phase,threaded,vol0,f,USEP)    
    end
end

function eos_g(wrapper::PTFlashWrapper{<:GammaPhi},p,T,z)
    RT = Rgas(wrapper)*T
    g_E = excess_gibbs_free_energy(__γ_unwrap(wrapper),p,T,z) #excess gibbs
    dg = -tpd_delta_g_vapour(wrapper,p,T,z)*RT #difference between gas fugacity and liquid activity, summed over z
    g_ideal = eos_g(idealmodel(wrapper),p,T,z) #ideal gibbs energy
    g0 = reference_state_eval(wrapper,p,T,z) #reference gibbs energy (equal to reference helmholtz energy)
    return g0 + g_ideal + g_E + dg
end

lb_volume(model::PTFlashWrapper{<:GammaPhi},T,z) = lb_volume(fluid_model(model),T,z)

function PT_property(model::PTFlashWrapper{<:GammaPhi},p,T,z,phase,threaded,vol0,f::F,USEP::Val{UseP}) where {F,UseP}
    if length(model) == 1
        return PT_property(fluid_model(model),p,T,z,phase,threaded,vol0,f,USEP)
    end
  
    if phase == :stable || is_unknown(phase)
        new_phase = identify_phase(model,p,T,z,vol0 = vol0)
        return PT_property(model,p,T,z,new_phase,threaded,vol0,f,USEP)
    end

    #shortcut for one-component models:

    #=
    Vapour properties are calculated with the fluid model
    Liquid properties are calculated via eos_g(PTFlashWrapper,p,T,z)
    =#
    if is_vapour(phase)
        return PT_property(gas_model(model),p,T,z,phase,threaded,vol0,f,USEP)
    elseif is_liquid(phase)
        if __γ_unwrap(model) isa IdealLiquidSolution
            #liquid phase + no activity: just delegate to the liquid model, whatever that model may be
            #even for saturated liquid volumes, you can get some props
            return PT_property(liquid_model(model),p,T,z,phase,threaded,vol0,f,USEP)
        end

        return PT_property_gibbs(model,p,T,z,f)
        #return PT_property(ActivityModelAresWrapper(model),p,T,z,phase,threaded,vol0,f,USEP)
    else
        throw(error("invalid phase specifier: $phase"))
    end
end

export GammaPhi
