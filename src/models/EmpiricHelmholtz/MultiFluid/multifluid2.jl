struct EmpiricMultiFluidParam <: EoSParam
    Mw::SingleParam{Float64}
    Tc::SingleParam{Float64}
    Pc::SingleParam{Float64}
    Vc::SingleParam{Float64}
    acentricfactor::SingleParam{Float64}
    lb_volume::SingleParam{Float64}
end

function EmpiricMultiFluidParam(components,pures)
    Mw = SingleParam("Mw",components,[pure.properties.Mw for pure in pures])
    Tc = SingleParam("Tc",components,[pure.properties.Tc for pure in pures])
    Pc = SingleParam("Pc",components,[pure.properties.Pc for pure in pures])
    Vc = SingleParam("Vc",components,1 ./ [pure.properties.rhoc for pure in pures])
    acentricfactor = SingleParam("acentric factor",components,[pure.properties.acentricfactor for pure in pures])
    lb_volume = SingleParam("lower bound volume",components,[pure.properties.lb_volume for pure in pures])
    return EmpiricMultiFluidParam(Mw,Tc,Pc,Vc,acentricfactor,lb_volume)
end

struct EmpiricMultiFluid{𝔸,𝕄} <: EmpiricHelmholtzModel
    components::Vector{String}
    params::EmpiricMultiFluidProperties
    pures::Vector{EmpiricSingleFluid{𝔸}}
    mixing::𝕄
    departure::ℙ
    Rgas::Float64
    references::Vector{String}
end

Rgas(model) = model.Rgas

function EmpiricMultiFluid(components;
    pure_userlocations = String[],
    mixing = EmpiricDeparture,
    mixing_userlocations = String[],
    estimate_pures = false,
    estimate_mixture = false,
    Rgas = nothing,
    verbose = false,
    )
    
    pures = [SingleFluid(comp;userlocations = pure_userlocations,verbose = verbose, xiang_deiters = estimate_pures, Rgas = Rgas) for comp in components]
    mixing = init_model(components,mixing,mixing_userlocations,verbose)
    departure = init_model(components,departure,departure_userlocations,verbose)
    
    #set Rgas
    if Rgas === nothing
    Rgases = [pure.properties.Rgas for pure in pures]
        R1 = first(Rgases)
        for Ri in Rgases
            if Ri != R1
                throw(error("unequal Rgas defined."))
            end
        end
        Rgas = R1
    else
        for i in eachindex(pures)
            pures[i] = __set_Rgas(pure[i],Rgas)
        end
    end

    params = EmpiricMultiFluidParam(components,pures)
    calculate_missing_mixing!(params,mixing)
    references = unique!(reduce(vcat,pure.references for pure in pures))
    return EmpiricMultiFluid(components,params,pures,mixing,departure,Rgas,references)
end

function reduced_delta(model,V,T,z,Σz = sum(z))
    Vᵣ = v_scale(model,V,T,z,∑z)
    Σz * Vᵣ/V
end

function reduced_tau(model,V,T,z,Σz = sum(z))
    Tᵣ = T_scale(model,V,T,z,∑z)
    Tᵣ / T
end

function a_ideal(model::EmpiricMultiFluid,V,T,z,∑z = sum(z))
    #log(δi) = log(ρ * vc[i]) = -log(V) + log(sum(z)) + log(vc[i])
    res = zero(V+T+first(z))
    m₀ = model.pures
    Tinv = 1/T
    Tc = model.params.Tc
    vc = model.params.Vc
    for i in 1:length(model)
        m₀ᵢ = m₀[i]
        a₀ᵢ = __get_k_alpha0(m₀ᵢ)*reduced_a_ideal(m₀ᵢ,Tc[i] * Tinv)
        zᵢ = z[i]
        res += zᵢ*a₀ᵢ
        res += xlogx(zᵢ,vc[i])
        res 
    end
    res /= ∑z
    res -= log(V)
    return res
end

function a_res(model::EmpiricMultiFluid,V,T,z)
    ∑z = sum(z)
    δ = reduced_delta(model,V,T,z,∑z)
    τ = reduced_tau(model,V,T,z,∑z)
    return multiparameter_a_res(model,V,T,z,model.departure,δ,τ,∑z)
end

function eos(model::EmpiricMultiFluid,V,T,z)
    ∑z = sum(z)
    a₀ = a_ideal(model,V,T,z,∑z)
    δ = v_scale(model,V,T,z,∑z)
    τ = T_scale(model,V,T,z,∑z)
    aᵣ = multiparameter_a_res(model,V,T,z,model.departure,δ,τ,∑z)
    return ∑z*@R̄()*T*(a₀+aᵣ)
end

function eos_res(model::EmpiricMultiFluid,V,T,z)
    ∑z = sum(z)
    δ = v_scale(model,V,T,z,∑z)
    τ = T_scale(model,V,T,z,∑z)
    aᵣ = multiparameter_a_res(model,V,T,z,model.departure,δ,τ,∑z)
    return ∑z*@R̄()*T*aᵣ
end

v_scale(model::EmpiricMultiFluid,V,T,z,∑z = sum(z)) = v_scale(model,V,T,z,model.mixing,∑z)
T_scale(model::EmpiricMultiFluid,V,T,z,∑z = sum(z)) = T_scale(model,V,T,z,model.mixing,∑z)

p_scale(model::EmpiricMultiFluid,z=SA[1.]) = dot(z,model.params.Pc.values)/sum(z)
 
T_scales(model::EmpiricMultiFluid,z=SA[1.]) = model.properties.Tc.values

#single functions, dispatch to pure
function x0_sat_pure(model::EmpiricMultiFluid,T)
    single_component_check(x0_sat_pure,model)
    x0_sat_pure(only(model.pures),T)
end

function x0_psat(model::EmpiricMultiFluid,T,crit = nothing)
    single_component_check(x0_psat,model)
    x0_psat(only(model.pures),T,crit)
end

function x0_saturation_temperature(model::EmpiricMultiFluid,p)
    single_component_check(x0_saturation_temperature,model)
    x0_saturation_temperature(only(model.pures),p)
end

function crit_pure(model::EmpiricMultiFluid)
    single_component_check(crit_pure,model)
    crit_pure(only(model.pures))
end

function lb_volume(model::EmpiricMultiFluid,z=SA[1.])
    return dot(z,model.params.lb_volume.values)
end

#use ideal gas
function x0_volume_gas(model::EmpiricMultiFluid,p,T,z=SA[1.])
    V = sum(z)*R̄*T/p
    return V
end

#use each available pure x0_volume_liquid
function x0_volume_liquid(model::EmpiricMultiFluid,T,z)
    v0 = zero(T+first(z))
    for (i,pure) in pairs(model.pures)
        v0 += z[i]*x0_volume_liquid(pure,T,SA[1.0])
    end
    return v0
end

function wilson_k_values(model::EmpiricMultiFluid,p,T,crit = nothing)
    n = length(model)
    K0 = zeros(typeof(p+T),n)
    pure = split_model.(model)
    _Tc = model.params.Tc.values
    _Pc = model.params.pc.values
    for i ∈ 1:n
        pure_i = pure[i]
        Tc,pc = _Tc[i],_Pc[i]
        ps = first(saturation_pressure(pure_i,0.7*Tc))
        ω = -log10(ps/pc) - 1.0
        K0[i] = exp(log(pc/p)+5.373*(1+ω)*(1-Tc/T))
    end
    return K0
end
