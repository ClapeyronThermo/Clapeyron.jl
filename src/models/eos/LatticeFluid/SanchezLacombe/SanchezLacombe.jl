struct SanchezLacombeParam <: EoSParam
    Mw::SingleParam{Float64}
    segment::SingleParam{Float64}
    epsilon::PairParam{Float64}
    vol::PairParam{Float64}
end



abstract type SanchezLacombeModel <: LatticeFluidModel end
include("mixing/mixing.jl")
struct SanchezLacombe{T <: SLMixingRule,I<:IdealModel} <:PRModel
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
    v_r,ε_r = mix_vε(model,V,T,z,mixing,r̄,Σz)
    r̄ = r̄/Σz
    v = V/Σz
    ρ̃ = r̄*v_r/v
    T̃ = R̄*T/ε_r
    _1 = one(V+T+first(z))
    return r̄*(-ρ̃ /T̃ + (_1/ρ̃  - _1)*log(1-ρ̃ )+_1)
end




