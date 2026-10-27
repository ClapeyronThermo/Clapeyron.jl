struct ChainLJScalingParam{T} <: ParametricEoSParam{T}
    Mw::SingleParam{T}
    segment::SingleParam{T}
    sigma::PairParam{T}
    epsilon::PairParam{T}
end

function ChainLJScalingParam(Mw,segment,sigma,epsilon)
    return build_parametric_param(ChainLJScalingParam,Mw,segment,sigma,epsilon)
end

@newmodelsimple ChainLJScaling LJScalingModel ChainLJScalingParam{T}

function ChainLJScaling(params::ChainLJScalingParam{T}) where T
    comps = params.Mw.components
    return ChainLJScaling(comps,params,comps) #use already allocated vector.
end

ChainLJScaling(param::EoSParam) = ChainLJScaling(ChainLJScalingParam(param.Mw,param.segment,param.sigma,param.epsilon))

function shape_factors(scaling::ChainLJScaling,ref,T,z)
    ϵ, σ, m = scaling.params.epsilon.values, scaling.params.sigma.values, scaling.params.segment.values
    σ3, ϵ̄, m̄ = σϵ_m_vdw1f(ϵ,σ,m,z)
    f = ϵ̄
    h = m̄*N_A*σ3
    return f,h
end

inverse_shape_factors(scaling::ChainLJScaling,ref,T,z) = shape_factors(scaling,ref,T,z)


function transform_params(::Type{ChainLJScaling},params,components)
    sigma = params["sigma"]
    sigma.values .*= 1E-10
    nc = length(components)
    segment = get!(params,"segment") do
        SingleParam("segment",components,ones(nc),fill(true,nc))
    end

    Mw = get!(params,"Mw") do
        SingleParam("Mw",components,zeros(nc),fill(true,nc))
    end    
    return saft_lorentz_berthelot(params)
end

"""
    LJScaling <: LJScalingModel
    LJScaling(components;
    userlocations = String[],
    verbose = false)

## Input parameters
- `sigma`: Single Parameter (`Float64`) - particle size `[Å]`
- `epsilon`: Single Parameter (`Float64`) - dispersion energy `[K]`
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Parameter (no units)
- `l`: Pair Parameter (`Float64`) (optional) - Binary Interaction Parameter (no units)
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)

## Model Parameters

- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `sigma`: Pair Parameter (`Float64`) - particle size `[m]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `epsilon`: Pair Parameter (`Float64`) - dispersion energy `[K]`

## Description

Lennard-Jones corresponding-states scaling model, used together with a Lennard-Jones reference equation of state.
The volume shape factor is additionally scaled by the weighted segment parameter
Their shape factors are defined as:
```
m̄    = Σᵢ xᵢ mᵢ
σ³   = Σᵢ Σⱼ xᵢ xⱼ σᵢⱼ³
ϵmix = (1/σ³) Σᵢ Σⱼ xᵢ xⱼ σᵢⱼ³ ϵᵢⱼ
σᵢⱼ  = (1 - lᵢⱼ)(σᵢ + σⱼ)/2        (Lorentz)
ϵᵢⱼ  = (1 - kᵢⱼ)√(ϵᵢ ϵⱼ)           (Berthelot)
f    = ϵmix
h    = m̄ N_A σ³
```
"""
ChainLJScaling

export ChainLJScaling