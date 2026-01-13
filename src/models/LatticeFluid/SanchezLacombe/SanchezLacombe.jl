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
    mixing::T
    params::SanchezLacombeParam
    idealmodel::I
    references::Array{String,1}
end

"""
    SanchezLacombe(components;
    idealmodel = BasicIdeal,
    mixing = SLk0k1lMixingRule,
    userlocations = String[],
    ideal_userlocations = String[],
    mixing_userlocations = String[],
    reference_state = false,
    verbose = false)

## Input parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `epsilon`: Single Parameter (`Float64`) - Nonbonded interaction energy per monomer `[J·mol⁻¹]`
- `vol`: Single Parameter (`Float64`) - Closed Packed Specific volume `[m³·mol⁻¹]`

## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `epsilon`: Pair Parameter (`Float64`) - Nonbonded interaction energy per monomer `[J·mol⁻¹]`
- `vol`: Pair Parameter (`Float64`) - Closed Packed Specific volume `[m³·mol⁻¹]`

## Input models
- `idealmodel`: Ideal Model
- `mixing`: Mixing model

## Description
Sanchez-Lacombe Lattice Fluid Equation of State.
```
xᵢ = zᵢ/∑zᵢ
r̄ = ∑xᵢrᵢ
vᵣ,εᵣ = mix_vε(model,V,T,z,model.mixing,r̄,∑zᵢ)
ρ̃ = r̄*vᵣ/v
T̃ = R̄*T/εᵣ
aᵣ = r̄*(- ρ̃ /T̃ + (1/ρ̃  - 1)*log(1 - ρ̃ ) + 1)
```

## References
1. Neau, E. (2002). A consistent method for phase equilibrium calculation using the Sanchez–Lacombe lattice–fluid equation-of-state. Fluid Phase Equilibria, 203(1–2), 133–140. [doi:10.1016/s0378-3812(02)00176-0](https://doi.org/10.1016/s0378-3812(02)00176-0)
"""
SanchezLacombe

const SL = SanchezLacombe

function SanchezLacombe(components;
    idealmodel = BasicIdeal,
    mixing = SLk0k1lMixingRule,
    userlocations = String[],
    ideal_userlocations = String[],
    mixing_userlocations = String[],
    reference_state = nothing,
    verbose = false)

    _components = format_components(components)
    params = getparams(_components, ["LatticeFluid/SanchezLacombe","properties/molarmass.csv"]; userlocations = userlocations, verbose = verbose)

    segment = params["segment"]
    unmixed_epsilon = params["epsilon"]
    unmixed_vol = params["vol"]
    Mw = params["Mw"]
    mixmodel = init_slmixing(mixing,_components,params,mixing_userlocations,verbose)
    ideal = init_model(idealmodel,_components,ideal_userlocations,verbose)
    premixed_vol,premixed_epsilon = sl_mix(unmixed_vol,unmixed_epsilon,mixmodel)
    packagedparams = SanchezLacombeParam(Mw, segment, premixed_epsilon, premixed_vol)
    references = ["10.1016/S0378-3812(02)00176-0"]
    model = SanchezLacombe(_components,mixmodel,packagedparams,ideal,references)
    set_reference_state!(model,reference_state;verbose)
    return model
end

function get_k(model::SanchezLacombe)
    return __SL_get_k(model,model.mixing)
end

function recombine_impl!(model::SanchezLacombe)
    recombine!(model.mixing)
    sl_mix!(model.params.vol,model.params.epsilon,model.mixing)
    return model
end

function get_l(model::SanchezLacombe)
    return __SL_get_l(model,model.mixing)
end

include("mixing/SLk0k1lrule.jl")
include("mixing/SLKrule.jl")

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
    #ρ̃/T̃ = ε_r*r̄*v_r/vRT
    #1/ρ̃ = v/(r̄*v_r)

    _1 = one(V+T+first(z))
    return r̄*(-ρ̃/T̃ + (_1/ρ̃  - _1)*log1p(-ρ̃)+_1)
end

function rmix(model::SanchezLacombe,V,T,z)
    r = model.params.segment.values
    r̄ = dot(z,r)/Σz
    return r̄
end

function lb_volume(model::SanchezLacombe,T,z)
    r = model.params.segment.values
    v = model.params.vol.values
    #v_r,ε_r = mix_vε(model,0.0,0.0,z,model.mixing,r̄,Σz)
    return sum(r[i]*z[i]*v[i,i] for i in @comps)
end

function T_scale(model::SanchezLacombe,z)
    Σz = sum(z)
    r = model.params.segment.values
    r̄ = dot(z,r)
    v_r,ε_r = mix_vε(model,0.0,0.0,z,model.mixing,r̄,Σz)
    return ε_r/R̄
end

function p_scale(model::SanchezLacombe,z)
    Σz = sum(z)
    r = model.params.segment.values
    r̄ = dot(z,r)
    v_r,ε_r = mix_vε(model,0.0,0.0,z,model.mixing,r̄,Σz)
    Ts = ε_r/R̄
    vs = v_r*r̄
    return R̄*Ts/vs
end

function x0_volume_liquid(model::SanchezLacombe,T,z)
    v_lb = lb_volume(model,T,z)
    return v_lb*1.1
end
#SL does not work with the virial coefficient
function x0_volume_gas(model::SanchezLacombe,p,T,z)
    return sum(z)*R̄*T/p
end
#=
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
    Tstar = Tr*369.89
    rhov = _propaneref_rhovsat(Tstar)
    vv = 1/rhov
    psat = pressure(model,vv,T)
    #vv = R̄*T/Ps
    vl = volume(model,psat,T,phase =:l)
    if isnan(vl)
        vv = nan
    end
    return (log10(vl),log10(vv))
end
=#
export SL,SanchezLacombe
