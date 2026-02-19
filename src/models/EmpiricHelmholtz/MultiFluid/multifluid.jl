struct MultiFluidParam <: EoSParam
    Mw::SingleParam{Float64}
    Tc::SingleParam{Float64}
    Pc::SingleParam{Float64}
    Vc::SingleParam{Float64}
    Tr::SingleParam{Float64}
    Vr::SingleParam{Float64}
    acentricfactor::SingleParam{Float64}
    lb_volume::SingleParam{Float64}
    reference_state::ReferenceState
end

function MultiFluidParam(components,pures,reference_state = nothing)
    Mw = SingleParam("Mw",components,[pure.properties.Mw for pure in pures])
    Tc = SingleParam("Tc",components,[pure.properties.Tc for pure in pures])
    Pc = SingleParam("Pc",components,[pure.properties.Pc for pure in pures])
    Vc = SingleParam("Vc",components,1 ./ [pure.properties.rhoc for pure in pures])
    Tr = SingleParam("Pc",components,[pure.properties.Tr for pure in pures])
    Vr = SingleParam("Vc",components,1 ./ [pure.properties.rhor for pure in pures])
    acentricfactor = SingleParam("acentric factor",components,[pure.properties.acentricfactor for pure in pures])
    lb_volume = SingleParam("lower bound volume",components,[pure.properties.lb_volume for pure in pures])
    ref = __init_reference_state_kw(reference_state)
    return MultiFluidParam(Mw,Tc,Pc,Vc,Tr,Vr,acentricfactor,lb_volume,ref)
end

struct MultiFluid{𝔸,𝕄,ℙ} <: EmpiricHelmholtzModel
    components::Vector{String}
    params::MultiFluidParam
    pures::Vector{SingleFluid{𝔸}}
    mixing::𝕄
    departure::ℙ
    Rgas::Float64
    references::Vector{String}
end

Rgas(model::MultiFluid) = model.Rgas

"""
    MultiFluid(components;
        idealmodel = nothing,
        ideal_userlocations = String[],
        pure_userlocations = String[],
        mixing = AsymmetricMixing,
        departure = EmpiricDeparture,
        mixing_userlocations = String[],
        departure_userlocations = String[],
        estimate_pure = false,
        estimate_mixing = :off,
        coolprop_userlocations = true,
        Rgas = nothing,
        reference_state = nothing,
         verbose = false)

## Input parameters
- JSON data (CoolProp and teqp format)

## Input models
- `idealmodel`: Ideal Model. If it is `nothing`, then it will parse the ideal model from the input JSON.
- `mixing`: mixing model for temperature and volume.
- `departure`: departure model

## Description

Instantiates a multi-component Empiric EoS model. `Rgas` can be used to set the value of the gas constant that is used during property calculations.

If `coolprop_userlocations` is true, then Clapeyron will try to look if the fluid is present in the CoolProp library.

If `estimate_pure` is true, then, if a JSON is not found, the pure model will be estimated, using the `XiangDeiters` model

`estimate_mixing` is used to fill missing mixing values in the case of using `AsymmetricMixing`. on other mixing models it has no effect.
 -  `estimate_mixing = :off` will perform no calculation of mixing parameter, throwing an error if missing values are found.
 -  `estimate_mixing = :lb` will perform Lorentz-Berthelot estimation of missing mixing parameters. (γT = βT = γv = βv = 1.0). additionally, you can pass `LorentzBerthelotMixing` to use `k` and `l` BIP instead.
 -  `estimate_mixing = :linear` will perform averaging of γT and γv so that `T(x) = ∑xᵢTᵢ` and `V(x) = ∑xᵢVᵢ` on missing mixing parameters. Additionally, you can use `LinearMixing` to perform this directly.

`Rgas` sets the value of the gas constant to be used by the multifluid. The default is the following:
- If `Rgas` is not specified and the input is a single component model, then the value of `Rgas` will be taken from the fluid json file.
- If `Rgas` is not specified and the input is a multi-component model, then the value of `Rgas` will be set to `Clapeyron.R̄ = Rgas() = 8.31446261815324` (2019 defined constant value)
"""
MultiFluid

function MultiFluid(components;
    idealmodel = nothing,
    pure_userlocations = String[],
    ideal_userlocations = String[],
    mixing = AsymmetricMixing,
    departure = EmpiricDeparture,
    mixing_userlocations = String[],
    departure_userlocations = String[],
    estimate_pure = false,
    estimate_mixing = :off,
    coolprop_userlocations = true,
    Rgas = nothing,
    reference_state = nothing,
    verbose = false)

    _components = format_components(components)
    idealmodels = if idealmodel === nothing
        fill(nothing,length(_components))
    else
        init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose,reference_state)
        split_model(init_idealmodel,1:length(_components))
    end

    pures = [
        SingleFluid(comp;
        userlocations = pure_userlocations,
        idealmodel = idealmodels[i],
        verbose = verbose,
        estimate_pure = estimate_pure,
        coolprop_userlocations = coolprop_userlocations,
        Rgas = Rgas
        )
        for (i,comp) in pairs(_components)]
    mixing = init_model(mixing,components,mixing_userlocations,verbose)
    departure = init_model(departure,components,departure_userlocations,verbose)
    params = MultiFluidParam(_components,pures,reference_state)
    references = unique!(reduce(vcat,pure.references for pure in pures))
    
    
    _Rgas = if Rgas == nothing
        if length(pures) != 1
            Clapeyron.Rgas()
        else
            Clapeyron.Rgas(pures[1])
        end
    else
        Rgas
    end
    model = MultiFluid(_components,params,pures,mixing,departure,_Rgas,references)
    recombine_mixing_reduced!(model,model.mixing,estimate_mixing)
    recombine_departure!(model,model.departure)
    set_reference_state!(model,verbose = verbose)
    return model
end

vT_scale(model::MultiFluid,V,T,z,Σz = sum(z)) = vT_scale(model,V,T,z,model.mixing,Σz)

function vT_scale(model,V,T,z,mixing,Σz)
    Vᵣ = v_scale(model,z,mixing,Σz)
    Tᵣ = T_scale(model,z,mixing,Σz)
    Vᵣ,Tᵣ
end

function reduced_delta_tau(model,V,T,z,Σz = sum(z))
    Vᵣ,Tᵣ = vT_scale(model,V,T,z,Σz)
    return Σz*Vᵣ/V, Tᵣ/T
end

function a_ideal(model::MultiFluid,V,T,z,∑z = sum(z))
    #log(δi) = log(ρ * vc[i]) = -log(V) + log(sum(z)) + log(vc[i])
    res = zero(V+T+first(z))
    m₀ = model.pures
    Tinv = 1/T
    Tc = model.params.Tr
    vc = model.params.Vr
    Rinv = 1/Rgas(model)
    for i in 1:length(model)
        m₀ᵢ = m₀[i]
        a₀ᵢ = reduced_a_ideal(m₀ᵢ,Tc[i] * Tinv)
        R₀ = m₀ᵢ.ideal.R0
        if !iszero(R₀)
            a₀ᵢ *=Rinv*R₀
        end
        zᵢ = z[i]
        res += zᵢ*a₀ᵢ
        res += xlogx(zᵢ,vc[i])
        res
    end
    res /= ∑z
    res -= log(V)
    return res
end

function a_res(model::MultiFluid,V,T,z)
    ∑z = sum(z)
    δ,τ = reduced_delta_tau(model,V,T,z,∑z)
    return multiparameter_a_res(model,V,T,z,model.departure,δ,τ,∑z)
end

function eos_impl(model::MultiFluid,V,T,z)
    ∑z = sum(z)
    a₀ = a_ideal(model,V,T,z,∑z)
    δ,τ = reduced_delta_tau(model,V,T,z,∑z)
    aᵣ = multiparameter_a_res(model,V,T,z,model.departure,δ,τ,∑z)
    return ∑z*Rgas(model)*T*(a₀+aᵣ) + reference_state_eval(model,V,T,z)
end

function eos_res(model::MultiFluid,V,T,z = SA[1.0])
    ∑z = sum(z)
    δ,τ = reduced_delta_tau(model,V,T,z,∑z)
    aᵣ = multiparameter_a_res(model,V,T,z,model.departure,δ,τ,∑z)
    return ∑z*Rgas(model)*T*aᵣ
end

v_scale(model::MultiFluid,z) = v_scale(model,z,sum(z))
T_scale(model::MultiFluid,z) = T_scale(model,z,sum(z))
v_scale(model::MultiFluid,z,∑z) = v_scale(model,z,model.mixing,∑z)
T_scale(model::MultiFluid,z,∑z) = T_scale(model,z,model.mixing,∑z)

p_scale(model::MultiFluid,z) = dot(z,model.params.Pc.values)/sum(z)

T_scales(model::MultiFluid,z=SA[1.]) = model.params.Tc.values

#single functions, dispatch to pure

function saturation_model(model::MultiFluid)
    return only(model.pures)
end

function lb_volume(model::MultiFluid,T,z)
    return dot(z,model.params.lb_volume.values)
end

function x0_crit_pure(model::MultiFluid,z)
    return (1.0,log10(v_scale(model,z)))
end


#use ideal gas
#function x0_volume_gas(model::MultiFluid,p,T,z)
#    
#end

has_fast_crit_pure(model::MultiFluid) = true

#use each available pure x0_volume_liquid
function x0_volume_liquid(model::MultiFluid,p,T,z)
    v0 = zero(Base.promote_eltype(model,p,T,z))
    lb_v = lb_volume(model,T,z)
    for (i,pure) in pairs(model.pures)
        if T > pure.properties.Tc
            lb_i = lb_volume(pure,T,SA[1.0])
            v0 += z[i]*1.01*lb_volume(pure,T,SA[1.0])
        else
            if p > x0_psat(pure,T)
                v0 += z[i]*1.01*lb_volume(pure,T,SA[1.0])
            else
                v0 += z[i]*x0_volume_liquid(pure,p,T,SA[1.0])
            end
        end
    end
    p0 = pressure(model,v0,T,z)
    for i in 1:10
        p0 > 0 && break
        v0 = 0.5v0 + 0.5*1.01*lb_v
        p0 = pressure(model,v0,T,z)
    end
    if p0 >= p
        return v0
    else
        return volume_bracket_refine(model,p,T,z,v0,lb_volume(model,T,z))
    end
end

function wilson_k_values!(K,model::MultiFluid,p,T,crit)
    n = length(model)
    pure = model.pures
    _Tc = model.params.Tc.values
    _Pc = model.params.Pc.values
    for i ∈ 1:n
        pure_i = pure[i]
        Tc,pc = _Tc[i],_Pc[i]
        ω = acentric_factor(pure_i,crit = (Tc,pc,NaN))
        K[i] = exp(log(pc/p)+ 5.3726985503194395*(1+ω)*(1-Tc/T))  #5.37 = log(10)*7/3
    end
    return K
end

function tp_flash_fast_K0!(K,model::MultiFluid,p,T,z)
    n = length(model)
    pure = model.pures
    _Tc = model.params.Tc.values
    _Pc = model.params.Pc.values
    for i ∈ 1:n
        pure_i = pure[i]
        Tc,pc = _Tc[i],_Pc[i]
        if T < Tc
            ps = x0_psat(pure_i,T)
            K[i] = ps/p
        else
            ps = x0_psat(pure_i,0.7*Tc)
            ω = -log10(ps/pc) - 1.0
            K[i] = exp(log(pc/p)+ 5.3726985503194395*(1+ω)*(1-Tc/T))
        end
    end
    return true
end

function split_pure_model(model::MultiFluid,splitter)
    pure_splitter = only.(splitter)
    model.pures[pure_splitter]
end

split_pure_model(model::MultiFluid,splitter::Int) = [model.pures[splitter]]
split_pure_model(model::MultiFluid,splitter::AbstractVector{<:Integer}) = model.pures[splitter]


#set reference states:
reference_state(model::MultiFluid) = model.params.reference_state
set_reference_state!(model::MultiFluid;verbose = false) = set_reference_state_empiric!(model;verbose)

function set_reference_state_empiric!(model;verbose = false)
    #handle cases where we don't need to do anything
    ref = reference_state(model)
    ref === nothing && return nothing
    ref.std_type == :no_set && return nothing
    if verbose
        @info "Calculating reference states for $model..."
        @info "Reference state type: $(info_color(ref.std_type))"
    end

    #allocate the appropiate caches.
    initialize_reference_state!(model,ref)
    pures = model.pures
    if all(iszero,ref.z0) #pure case
        pure_refs = split_model(ref,1:length(model))
        _set_reference_state!.(pures,SA[1.0],pure_refs)
        ref.a0 .= only.(getfield.(pure_refs,:a0))
        ref.a1 .= only.(getfield.(pure_refs,:a1))
    else
        _set_reference_state!(model,ref.z0)
    end
    for (i,pure) in pairs(pures)
        ref_a = pure.ideal.ref_a
        ref_a[1] = pure_refs[i].a0[1]
        ref_a[2] = pure_refs[i].a1[1]
    end
    return model
end

export MultiFluid

