struct PeTSParam <: EoSParam
    Mw::SingleParam{Float64}
    segment::SingleParam{Float64}
    sigma::PairParam{Float64}
    epsilon::PairParam{Float64}
end
abstract type PeTSModel <: EoSModel end

struct PeTS{T <: IdealModel} <:PeTSModel
    components::Array{String,1}
    params::PeTSParam
    idealmodel::T
    references::Array{String,1}
end

@registermodel PeTS
"""
    PeTSModel <: EoSModel

    PeTS(components; 
    idealmodel=BasicIdeal,
    userlocations=String[],
    ideal_userlocations=String[],
    verbose=false)

## Input parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Single Parameter (`Float64`) - Segment Diameter [`A°`]
- `epsilon`: Single Parameter (`Float64`) - Reduced dispersion energy  `[K]`
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Paramater (no units)

## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Pair Parameter (`Float64`) - Mixed segment Diameter `[m]`
- `epsilon`: Pair Parameter (`Float64`) - Mixed reduced dispersion energy`[K]`

## Input models
- `idealmodel`: Ideal Model

## Description

Perturbed, Truncated and Shifted (PeTS) Equation of State.

## References
1. Heier, M., Stephan, S., Liu, J., Chapman, W. G., Hasse, H., & Langenbach, K. (2018). Equation of state for the Lennard-Jones truncated and shifted fluid with a cut-off radius of 2.5 σ based on perturbation theory and its applications to interfacial thermodynamics. Molecular Physics, 116(15–16), 2083–2094. [doi:10.1080/00268976.2018.1447153](https://doi.org/10.1080/00268976.2018.1447153)
"""
PeTS

export PeTS

function PeTS(components;
    idealmodel=BasicIdeal,
    userlocations=String[],
    ideal_userlocations=String[],
    verbose=false)
    params,sites = getparams(components, ["SAFT/PCSAFT","properties/molarmass.csv"]; userlocations=userlocations, verbose=verbose)
    
    segment = params["segment"]
    k = get(params,"k",nothing)
    Mw = params["Mw"]
    params["sigma"].values .*= 1E-10
    sigma = sigma_LorentzBerthelot(params["sigma"])
    epsilon = epsilon_LorentzBerthelot(params["epsilon"], k)
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)

    packagedparams = PeTSParam(Mw, segment, sigma, epsilon)
    references = ["10.1080/00268976.2018.1447153"]
    model = PeTS(components,packagedparams,init_idealmodel,references)
    return model
end

function a_res(model::PeTSModel, V, T, z)
    _data = @f(data)
    return @f(a_ref,_data) + @f(a_pert,_data) 
end

function data(model::PeTSModel,V,T,z)
    σ3,ϵ̄,m̄ = σϵ_m_vdw1f(model,V,T,z)
    T̃ = T/ϵ̄
    N  = N_A*sum(z)
    ρS = N/V*m̄
    ρ̃  = ρS*σ3
    d̃ = d_pets(T̃)
    η = π*ρ̃ *d̃^3 / 6
    return (η,ρ̃ ,T̃,m̄)
end


function d(model::PeTSModel, V, T, z,_data=@f(data))
    η,ρ̃ ,T̃ = _data
    return d_pets(T̃)
end

function d_pets(T̃)
    #return 1 - 0.127112544*exp(-3.052785558/T̃)
    #log(0.127112544) = -2.0626824117148774
    return -expm1(-3.052785558/T̃ - 2.0626824117148774)
end

a_hs(model::PeTSModel,V,T,z,_data=@f(data)) = a_ref(model,V,T,z,_data) 

function a_ref(model::PeTSModel, V, T, z,_data=@f(data))
    η,ρ̃ ,T̃,_ = _data
    return η*(4-3η)/(1-η)^2
end

function a_pert(model::PeTSModel, V, T, z,_data=@f(data))
    η,ρ̃ ,T̃,_ = _data
    I1 = evalpoly(η,PeTS_A)
    I2 = evalpoly(η,PeTS_B)
    ã1 = -2*π*ρ̃ *I1/ T̃
    ã2 = -π*ρ̃ *I2*(1 + 2*η*(4 - η)/(1 - η)^4)^-1 / T̃^2
    return (ã1 + ã2)
end

const PeTS_A = (
    0.690603404,
    1.189317012,
    1.265604153,
    -24.34554201,
    93.67300357,
    -157.8773415,
    96.93736697)

const PeTS_B = (
    0.664852128,
    2.10733079,
    -9.597951213,
    -17.37871193,
    30.17506222,
    209.3942909,
    -353.2743581)

function lb_volume(model::PeTSModel,z=SA[1.0])
    σ3,_,m̄ = σϵ_m_vdw1f(model,1.0,1.0,z)
    return m̄*N_A*σ3*π/6
end

function x0_volume_liquid(model::PeTSModel,T,z)
    v_lb = lb_volume(model,z)
    return v_lb*1.8
end

function p_scale(model::PeTSModel,z = SA[1.0])
    σ3,ϵ,m̄ = σϵ_m_vdw1f(model,1.0,1.0,z)
    v = m̄*N_A*σ3
    return R̄*ϵ/v
end

function T_scale(model::PeTSModel,z=SA[1.0])
    σ3,ϵ,m̄ = σϵ_m_vdw1f(model,1.0,1.0,z)
    return ϵ
end

function x0_crit_pure(model::PeTSModel)
    lb_v = lb_volume(model)
    (1.08, log10(lb_v/0.32))
end