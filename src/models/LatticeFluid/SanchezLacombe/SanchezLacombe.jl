struct SanchezLacombeParam <: EoSParam
    Mw::SingleParam{Float64}
    segment::SingleParam{Float64}
    epsilon::PairParam{Float64}
    vol::PairParam{Float64}
end



abstract type SanchezLacombeModel <: LatticeFluidModel end
include("mixing/mixing.jl")
struct SanchezLacombe{T <: SLMixingRule,I<:IdealModel} <:SanchezLacombeModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    mixing::T
    params::SanchezLacombeParam
    idealmodel::I
    references::Array{String,1}
end
@registermodel SanchezLacombe
const SL = SanchezLacombe


function SanchezLacombe(components; 
    idealmodel=BasicIdeal, 
    mixing = SLk0k1lMixingRule, 
    userlocations=String[], 
    ideal_userlocations=String[], 
    mixing_userlocations = String[],
    verbose=false)
    params = getparams(components, ["LatticeFluid/SanchezLacombe","properties/molarmass.csv"]; userlocations=userlocations, verbose=verbose)
    
    segment = params["segment"]
    unmixed_epsilon = params["epsilon"]
    unmixed_vol = params["vol"]
    unmixed_epsilon.values #.*= k_B #to convert from temperature to eps
    unmixed_vol.values .*= 1e-6 #convert from cm3/mol to m3/mol
    Mw = params["Mw"]
    mixmodel = init_model(mixing,components,mixing_userlocations,verbose)
    ideal = init_model(idealmodel,components,ideal_userlocations,verbose)
    premixed_vol,premixed_epsilon = sl_mix(unmixed_vol,unmixed_epsilon,mixmodel)
    packagedparams = SanchezLacombeParam(Mw, segment, premixed_epsilon, premixed_vol)
    references = ["10.1016/S0378-3812(02)00176-0"]
    icomponents = 1:length(components)
    model = SanchezLacombe(components,icomponents,mixmodel,packagedparams,ideal,references)
    return model
end

include("mixing/SLk0k1lrule.jl")

function a_res(model::SanchezLacombe,V,T,z=SA[1.0])
    Σz = sum(z)     
    r = model.params.segment.values
    mixing = model.mixing
    r̄ = dot(z,r)
    r̄ = r̄/Σz
    v_r,ε_r = mix_vε(model,V,T,z,mixing,r̄,Σz)
    v = V/Σz
    ρ̃ = r̄*v_r/v
    T̃ = R̄*T/ε_r
    _1 = one(V+T+first(z))
    return r̄*(-ρ̃ /T̃ + (_1/ρ̃  - _1)*log(1-ρ̃ )+_1)
end

function lb_volume(model::SanchezLacombe,z=SA[1.0])
    Σz = sum(z)
    r = model.params.segment.values
    v = model.params.vol.diagvalues
    r̄ = dot(z,r)
    #v_r,ε_r = mix_vε(model,0.0,0.0,z,model.mixing,r̄,Σz)
    return sum(r[i]*z[i]*v[i] for i in @comps)
end

function T_scale(model::SanchezLacombe,z=SA[1.0])
    Σz = sum(z)
    r = model.params.segment.values
    r̄ = dot(z,r)
    v_r,ε_r = mix_vε(model,0.0,0.0,z,model.mixing,r̄,Σz)
    return ε_r/R̄
end

function p_scale(model::SanchezLacombe,z=SA[1.0])
    Σz = sum(z)
    r = model.params.segment.values
    r̄ = dot(z,r)
    v_r,ε_r = mix_vε(model,0.0,0.0,z,model.mixing,r̄,Σz)
    Ts = ε_r/R̄
    vs = v_r*r̄
    return R̄*Ts/vs
end

function x0_volume_liquid(model::SanchezLacombe,T,z)
    v_lb = lb_volume(model,z)
    return v_lb*1.1
end
#SL does not work with the virial coefficient
function x0_volume_gas(model::SanchezLacombe,p,T,z)
    return sum(z)*R̄*T/p
end

function x0_sat_pure(model::SanchezLacombe,T,z=SA[1.0])
    Σz = sum(z)
    r = model.params.segment.values
    r̄ = dot(z,r)
    v_r,ε_r = mix_vε(model,0.0,T,z,model.mixing,r̄,Σz)
    Ts = ε_r/R̄
    vs = v_r*r̄
    Ps = R̄*Ts/vs
    Tr = T/Ts
    nan = zero(Tr)/zero(Tr)
    Tr > 1 && return [nan,nan]
    Tstar = Tr*PropaneRef_consts.T_c
    rhov = _propaneref_rhovsat(Tstar)
    vv = 1/rhov
    psat = pressure(model,vv,T)
    #vv = R̄*T/Ps
    vl = volume(model,psat,T,phase =:l)
    if isnan(vl)
        vv = nan
    end
    return [log10(vl),log10(vv)]
end

export SL,SanchezLacombe,SLk0k1lMixingRule
