struct MultiFluidParam <: EoSParam
    Mw::SingleParam{Float64}
    Tc::SingleParam{Float64}
    Pc::SingleParam{Float64}
    Vc::SingleParam{Float64}
    acentricfactor::SingleParam{Float64}
    lb_volume::SingleParam{Float64}
end

function MultiFluidParam(components,pures)
    Mw = SingleParam("Mw",components,[pure.properties.Mw for pure in pures])
    Tc = SingleParam("Tc",components,[pure.properties.Tc for pure in pures])
    Pc = SingleParam("Pc",components,[pure.properties.Pc for pure in pures])
    Vc = SingleParam("Vc",components,1 ./ [pure.properties.rhoc for pure in pures])
    acentricfactor = SingleParam("acentric factor",components,[pure.properties.acentricfactor for pure in pures])
    lb_volume = SingleParam("lower bound volume",components,[pure.properties.lb_volume for pure in pures])
    return MultiFluidParam(Mw,Tc,Pc,Vc,acentricfactor,lb_volume)
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
        Rgas = R̄,
        verbose = false,
        )
## Input parameters
- JSON data (CoolProp and teqp format)

## Input models
- `idealmodel`: Ideal Model. if it is `nothing`, then it will parse the ideal model from the input JSON.
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

`Rgas` sets the value of the gas constant to be used by the multifluid. defaults to `Clapeyron.R̄ = Rgas() = 8.31446261815324` (2019 defined constant value)
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
    Rgas = R̄,
    verbose = false,
    )

    _components = format_components(components)
    if idealmodel === nothing
        idealmodels = FillArrays.Fill(nothing,length(_components))
    else
        init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
        idealmodels = split_model(init_idealmodel,1:length(_components))
    end

    pures = [
        SingleFluid(comp;
        userlocations = pure_userlocations,
        idealmodel = idealmodels[i],
        verbose = verbose,
        estimate_pure = estimate_pure,
        coolprop_userlocations = coolprop_userlocations,
        )
        for (i,comp) in pairs(_components)]
    mixing = init_model(mixing,components,mixing_userlocations,verbose)
    departure = init_model(departure,components,departure_userlocations,verbose)
    params = MultiFluidParam(_components,pures)
    references = unique!(reduce(vcat,pure.references for pure in pures))
    model = MultiFluid(_components,params,pures,mixing,departure,Rgas,references)
    recombine_mixing!(model,model.mixing,estimate_mixing)
    recombine_departure!(model,model.departure)
    return model
end

function reduced_delta(model,V,T,z,Σz = sum(z))
    Vᵣ = v_scale(model,z,Σz)
    Σz * Vᵣ/V
end

function reduced_tau(model,V,T,z,Σz = sum(z))
    Tᵣ = T_scale(model,z,Σz)
    Tᵣ / T
end

function a_ideal(model::MultiFluid,V,T,z,∑z = sum(z))
    #log(δi) = log(ρ * vc[i]) = -log(V) + log(sum(z)) + log(vc[i])
    res = zero(V+T+first(z))
    m₀ = model.pures
    Tinv = 1/T
    Tc = model.params.Tc
    vc = model.params.Vc
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
    δ = reduced_delta(model,V,T,z,∑z)
    τ = reduced_tau(model,V,T,z,∑z)
    return multiparameter_a_res(model,V,T,z,model.departure,δ,τ,∑z)
end

function eos(model::MultiFluid,V,T,z = SA[1.0])
    ∑z = sum(z)
    a₀ = a_ideal(model,V,T,z,∑z)
    δ = reduced_delta(model,V,T,z,∑z)
    τ = reduced_tau(model,V,T,z,∑z)
    aᵣ = multiparameter_a_res(model,V,T,z,model.departure,δ,τ,∑z)
    return ∑z*@R̄()*T*(a₀+aᵣ)
end

function eos_res(model::MultiFluid,V,T,z = SA[1.0])
    ∑z = sum(z)
    δ = reduced_delta(model,V,T,z,∑z)
    τ = reduced_tau(model,V,T,z,∑z)
    aᵣ = multiparameter_a_res(model,V,T,z,model.departure,δ,τ,∑z)
    return ∑z*@R̄()*T*aᵣ
end

v_scale(model::MultiFluid,z = SA[1.0],∑z = sum(z)) = v_scale(model,z,model.mixing,∑z)
T_scale(model::MultiFluid,z = SA[1.0],∑z = sum(z)) = T_scale(model,z,model.mixing,∑z)

p_scale(model::MultiFluid,z=SA[1.]) = dot(z,model.params.Pc.values)/sum(z)

T_scales(model::MultiFluid,z=SA[1.]) = model.properties.Tc.values

#single functions, dispatch to pure
function x0_sat_pure(model::MultiFluid,T)
    single_component_check(x0_sat_pure,model)
    x0_sat_pure(only(model.pures),T)
end

function x0_psat(model::MultiFluid,T,crit = nothing)
    single_component_check(x0_psat,model)
    x0_psat(only(model.pures),T,crit)
end

function x0_saturation_temperature(model::MultiFluid,p)
    single_component_check(x0_saturation_temperature,model)
    x0_saturation_temperature(only(model.pures),p)
end

function crit_pure(model::MultiFluid)
    single_component_check(crit_pure,model)
    crit_pure(only(model.pures))
end

function lb_volume(model::MultiFluid,z=SA[1.])
    return dot(z,model.params.lb_volume.values)
end

#use ideal gas
function x0_volume_gas(model::MultiFluid,p,T,z=SA[1.])
    V = sum(z)*R̄*T/p
    return V
end

#use each available pure x0_volume_liquid
function x0_volume_liquid(model::MultiFluid,T,z)
    v0 = zero(T+first(z))
    for (i,pure) in pairs(model.pures)
        v0 += z[i]*x0_volume_liquid(pure,T,SA[1.0])
    end
    return v0
end

function wilson_k_values!(K,model::MultiFluid,p,T,crit = nothing)
    n = length(model)
    pure = split_model.(model)
    _Tc = model.params.Tc.values
    _Pc = model.params.Pc.values
    for i ∈ 1:n
        pure_i = pure[i]
        Tc,pc = _Tc[i],_Pc[i]
        ps = first(saturation_pressure(pure_i,0.7*Tc))
        ω = -log10(ps/pc) - 1.0
        K[i] = exp(log(pc/p)+5.373*(1+ω)*(1-Tc/T))
    end
    return K
end

export MultiFluid

