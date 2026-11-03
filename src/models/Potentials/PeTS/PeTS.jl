


abstract type PeTSModel <: LJCorrespondingStatesModel end

struct PeTS{T,I} <: PeTSModel
    components::Vector{String}
    params::ChainLJScalingParam{T}
    idealmodel::I
    references::Vector{String}
    molecular_units::Bool
end

default_locations(::Type{PeTS}) = ["SAFT/PCSAFT","properties/molarmass.csv"]
default_references(::Type{PeTS}) = ["10.1080/00268976.2018.1447153"]

"""
    PeTSModel <: EoSModel

    PeTS(components;
    idealmodel = BasicIdeal,
    userlocations = String[],
    ideal_userlocations = String[],
    reference_state = nothing,
    molecular_units = false,
    verbose = false)

## Input parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Single Parameter (`Float64`) - Segment Diameter `[Å]`
- `epsilon`: Single Parameter (`Float64`) - Reduced dispersion energy `[K]`
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Parameter (no units)

## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Pair Parameter (`Float64`) - Mixed segment Diameter `[m]`
- `epsilon`: Pair Parameter (`Float64`) - Mixed reduced dispersion energy`[K]`

## Input models
- `idealmodel`: Ideal Model

## Description

Perturbed, Truncated and Shifted (PeTS) Equation of State. Scaled via [`LJScaling`](@ref)

The `molecular_units` argument sets if the scaling is done in SI units or molecular units.
On SI units, the input sigma diameter is scaled from `[Å]` to `[m]` (once, at construction time), the volume is scaled by the avogadro number, and the gas constant is set to `R = 8.314...`;
on molecular units, no additional scaling is applied and `R = 1.0`

## References
1. Heier, M., Stephan, S., Liu, J., Chapman, W. G., Hasse, H., & Langenbach, K. (2018). Equation of state for the Lennard-Jones truncated and shifted fluid with a cut-off radius of 2.5 σ based on perturbation theory and its applications to interfacial thermodynamics. Molecular Physics, 116(15–16), 2083–2094. [doi:10.1080/00268976.2018.1447153](https://doi.org/10.1080/00268976.2018.1447153)
"""
function PeTS(components;
    idealmodel = BasicIdeal,
    userlocations = String[],
    ideal_userlocations = String[],
    molecular_units = false,
    reference_state = nothing,
    verbose = false)

    _components = format_components(components)
    params1 = getparams(_components, default_locations(PeTS); userlocations = userlocations, verbose = verbose)
    params2 = transform_params_chainljscaling(params1,_components,molecular_units)
    packagedparams = build_eosparam(ChainLJScalingParam,params2)
    T = eltype(packagedparams)
    init_idealmodel = init_model(idealmodel,_components,ideal_userlocations,verbose)
    references = default_references(PeTS)
    model = PeTS{T,typeof(init_idealmodel)}(_components, packagedparams, init_idealmodel, references, molecular_units)
    set_reference_state!(model,reference_state;verbose)
    return model
end

"""
    PeTSLJ <: LJPotentialModel
    PeTSLJ()

## Description

Perturbed, Truncated and Shifted (PeTS) Equation of State. Unscaled model.

## References
1. Heier, M., Stephan, S., Liu, J., Chapman, W. G., Hasse, H., & Langenbach, K. (2018). Equation of state for the Lennard-Jones truncated and shifted fluid with a cut-off radius of 2.5 σ based on perturbation theory and its applications to interfacial thermodynamics. Molecular Physics, 116(15–16), 2083–2094. [doi:10.1080/00268976.2018.1447153](https://doi.org/10.1080/00268976.2018.1447153)
"""
struct PeTSLJ <: LJPotentialModel end

reference_model(::PeTSModel) = PeTSLJ()
scaling_model(model::PeTSModel) = ChainLJScaling(model.params,model.molecular_units)
default_ignore_missing_singleparams(model::PeTSModel) =  ["segment"]

function a_potential(model::PeTSLJ,V,T)
    z = SA[1.0]
    _data = @f(data)
    return @f(a_ref,_data) + @f(a_pert,_data)
end

function data(model::PeTSLJ,V,T̃,z)
    d̃ = d_pets(T̃)
    ρ̃  = 1/V
    η = π*ρ̃ *d̃^3 / 6
    return (η,ρ̃ ,T̃)
end

function d(model::PeTSLJ, V, T, z,_data=@f(data))
    η,ρ̃ ,T̃ = _data
    return d_pets(T̃)
end

function d_pets(T̃)
    #return 1 - 0.127112544*exp(-3.052785558/T̃)
    #log(0.127112544) = -2.0626824117148774
    return -expm1(-3.052785558/T̃ - 2.0626824117148774)
end

a_hs(model::PeTSLJ,V,T,z,_data=@f(data)) = a_ref(model,V,T,z,_data)

function a_ref(model::PeTSLJ, V, T, z,_data=@f(data))
    η,ρ̃ ,T̃ = _data
    return η*(4-3η)/(1-η)^2
end

function a_pert(model::PeTSLJ, V, T, z,_data=@f(data))
    η,ρ̃ ,T̃ = _data
    I1 = evalpoly(η,PeTS_A)
    I2 = evalpoly(η,PeTS_B)
    ã1 = -2*π*ρ̃ *I1/ T̃
    ã2 = -π*ρ̃ *I2*(1 + 2*η*(4 - η)/(1 - η)^4)^-1 / T̃^2
    return (ã1 + ã2)
end

function lb_volume(model::PeTSLJ,T,z)
    d̃ = d_pets(T)
    return sum(z)*π*d̃^3 / 6
end

function x0_volume_liquid(model::PeTSLJ,p,T,z)
    return 1.8*lb_volume(model,T,z)
end

crit_pure(model::PeTSLJ) = (1.0890488908356557, 0.10204180775002869, 3.234497969253909)
has_fast_crit_pure(::PeTSLJ) = true

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

export PeTS, PeTSLJ