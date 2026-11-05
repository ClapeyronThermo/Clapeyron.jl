struct ChainLJScalingParam{T} <: ParametricEoSParam{T}
    Mw::SingleParam{T}
    segment::SingleParam{T}
    sigma::PairParam{T}
    epsilon::PairParam{T}
end

function ChainLJScalingParam(Mw,segment,sigma,epsilon)
    return build_parametric_param(ChainLJScalingParam,Mw,segment,sigma,epsilon)
end

struct ChainLJScaling{T} <: LJScalingModel
    components::Vector{String}
    params::ChainLJScalingParam{T}
    molecular_units::Bool
end

function ChainLJScaling(params::ChainLJScalingParam{T},molecular_units::Bool) where T
    comps = params.Mw.components
    return ChainLJScaling(comps,params,molecular_units) #use already allocated vector.
end

ChainLJScaling(param::EoSParam) = ChainLJScaling(param,false)
ChainLJScaling(param::EoSParam,molecular_units::Bool) = ChainLJScaling(ChainLJScalingParam(param.Mw,param.segment,param.sigma,param.epsilon),molecular_units)

function shape_factors(scaling::ChainLJScaling,ref,T,z)
    ϵ, σ, m = scaling.params.epsilon.values, scaling.params.sigma.values, scaling.params.segment.values
    σ3, ϵ̄, m̄ = σϵ_m_vdw1f(ϵ,σ,m,z)
    f = ϵ̄
    h = m̄*σ3
    if !scaling.molecular_units
        h *= N_A
    end
    return f,h
end

inverse_shape_factors(scaling::ChainLJScaling,ref,T,z) = shape_factors(scaling,ref,T,z)

function transform_params_chainljscaling(params,components,molecular_units)
    params2 = transform_params_ljscaling(params,components,molecular_units)

    get!(params2,"segment") do
        nc = length(components)
        SingleParam("segment",components,ones(nc),fill(true,nc))
    end

    return params2
end

"""
    LJScaling <: LJScalingModel
    LJScaling(components;
    userlocations = String[],
    molecular_units = false,
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
h    = m̄σ³ * scaling
```
The `molecular_units` argument sets if the scaling is done in SI units or molecular units.
On SI units, the input sigma diameter is scaled from `[Å]` to `[m]` (once, at construction time), the volume is scaled by the avogadro number, and the gas constant is set to `R = 8.314...`;
on molecular units, no additional scaling is applied and `R = 1.0`
"""
function ChainLJScaling(components;
                userlocations = String[],
                verbose = false,
                molecular_units = false)


    _components = format_components(components)
    params1 = getparams(_components, String[]; userlocations = userlocations, verbose = verbose)
    params2 = transform_params_chainljscaling(params1,_components,molecular_units)
    packagedparams = build_eosparam(ChainLJScalingParam,params2)
    T = eltype(packagedparams)
    model = ChainLJScaling{T}(_components, packagedparams, molecular_units)
    return model
end

Rgas(model::ChainLJScaling) = ifelse(model.molecular_units,1.0,R̄)

export ChainLJScaling