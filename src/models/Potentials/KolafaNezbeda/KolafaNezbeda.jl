
abstract type KolafaNezbedaModel <: LJCorrespondingStatesModel end

struct KolafaNezbeda{T,I} <: KolafaNezbedaModel
    components::Vector{String}
    params::LJScalingParam{T}
    idealmodel::I
    references::Vector{String}
    molecular_units::Bool
end
function show_info(io,model::KolafaNezbeda) 
    molecular_units = model.molecular_units
    if molecular_units
        println(io)
        print(io," molecular units: ",molecular_units)
    end
end

default_locations(::Type{KolafaNezbeda}) = ["Potentials/KolafaNezbeda"]
default_references(::Type{KolafaNezbeda}) = ["10.1016/0378-3812(94)02573-G"]

reference_model(::KolafaNezbedaModel) = KolafaNezbedaLJ()
scaling_model(model::KolafaNezbedaModel) = LJScaling(model.params,model.molecular_units)
default_ignore_missing_singleparams(model::KolafaNezbedaModel) =  ["segment"]

"""
    KolafaNezbedaModel <: EoSModel

    KolafaNezbeda(components;
    idealmodel = BasicIdeal,
    userlocations = String[],
    ideal_userlocations = String[],
    reference_state = nothing,
    molecular_units = false,
    verbose = false)

## Input parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `sigma`: Single Parameter (`Float64`) - Lennard-Jones size parameter `[Å]`
- `epsilon`: Single Parameter (`Float64`) - Lennard-Jones energy parameter `[K]`
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Parameter (no units)

## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `sigma`: Pair Parameter (`Float64`) - Mixed Lennard-Jones size parameter `[m]`
- `epsilon`: Pair Parameter (`Float64`) - Mixed Lennard-Jones energy parameter `[K]`

## Input models
- `idealmodel`: Ideal Model

## Description

Kolafa-Nezbeda equation of state for the Lennard-Jones fluid. Scaled via [`LJScaling`](@ref)

The `molecular_units` argument sets if the scaling is done in SI units or molecular units.
On SI units, the input sigma diameter is scaled from `[Å]` to `[m]` (once, at construction time), the volume is scaled by the avogadro number, and the gas constant is set to `R = 8.314...`;
on molecular units, no additional scaling is applied and `R = 1.0`

## References
1. Kolafa, J., & Nezbeda, I. (1994). The Lennard-Jones fluid: An accurate analytic and theoretically-based equation of state. Fluid Phase Equilibria, 100, 1–34. [doi:10.1016/0378-3812(94)02573-G](https://doi.org/10.1016/0378-3812(94)02573-G)
"""
function KolafaNezbeda(components;
    idealmodel = BasicIdeal,
    userlocations = String[],
    ideal_userlocations = String[],
    molecular_units = false,
    reference_state = nothing,
    verbose = false)

    _components = format_components(components)
    params1 = getparams(_components, default_locations(KolafaNezbeda); userlocations = userlocations, verbose = verbose)
    params2 = transform_params_ljscaling(params1,_components,molecular_units)
    packagedparams = build_eosparam(LJScalingParam,params2)
    T = eltype(packagedparams)
    init_idealmodel = init_model(idealmodel,_components,ideal_userlocations,verbose)
    references = default_references(KolafaNezbeda)
    model = KolafaNezbeda{T,typeof(init_idealmodel)}(_components, packagedparams, init_idealmodel, references, molecular_units)
    set_reference_state!(model,reference_state;verbose)
    return model
end

#nonlinear parameter γ of the damping function, eq. (15)
const KN_γ = 1.92907278
const KN_ΔB2 = (0.02459877, 0.0, -7.02181962, 2.90616279, -4.13749995, 0.87361369, 0.43102052, -0.58544978)

#residual term Cᵢⱼ, table 3. one tuple per i, coefficients ordered in j = 2:6
const KN_C0 = (2.01546797, -28.17881636, 28.28313847, -10.42402873)
const KN_C1 = (-19.58371655, 75.62340289, -120.70586598, 93.92740328, -27.37737354)
const KN_C2 = (29.34470520, -112.35356937, 170.64908980, -123.06669187, 34.42288969)
const KN_C4 = (-13.37031968, 65.38059570, -115.09233113, 88.91973082, -25.62099890)

@newmodelsingleton KolafaNezbedaLJ LJPotentialModel

function a_potential(model::KolafaNezbedaLJ, V, T)
    z = SA[1.0]
    _data = @f(data)
    return @f(a_hs,_data) + @f(a_pert,_data)
end

function data(model::KolafaNezbedaLJ,V,T,z)
    ρ̃ = sum(z)/V
    d̃ = d_kn(T)
    η = π*ρ̃*d̃^3/6
    return (η,ρ̃,T,d̃)
end

#hybrid Barker-Henderson hard sphere diameter
function d_kn(T̃)
    sT = sqrt(T̃)
    return 1.080142248 + 0.011117524/T̃ - 0.076383859/sT + 0.000693129*sT - 0.063920968*log(T̃)
end

#residual second virial coefficient
function ΔB2_kn(T̃)
    return evalpoly(1/sqrt(T̃),KN_ΔB2)
end

#hard sphere term, associated to the Boublik-Nezbeda hard sphere EoS
function a_hs(model::KolafaNezbedaLJ, V, T, z,_data = @f(data))
    η,ρ̃,T̃,d̃ = _data
    return 5/3*log1p(-η) + η*(34 - 33*η + 4*η*η)/(6*(1-η)^2)
end

function a_pert(model::KolafaNezbedaLJ, V, T, z,_data = @f(data))
    return @f(a_B2,_data) + @f(a_corr,_data)
end

function a_B2(model::KolafaNezbedaLJ, V, T, z,_data = @f(data))
    η,ρ̃,T̃,d̃ = _data
    return ρ̃*ΔB2_kn(T̃)*exp(-KN_γ*ρ̃*ρ̃)
end

function a_corr(model::KolafaNezbedaLJ, V, T, z,_data = @f(data))
    η,ρ̃,T̃,d̃ = _data
    u = 1/sqrt(T̃)
    C0 = evalpoly(ρ̃,KN_C0)
    C1 = evalpoly(ρ̃,KN_C1)
    C2 = evalpoly(ρ̃,KN_C2)
    C4 = evalpoly(ρ̃,KN_C4)
    return ρ̃*ρ̃*u*u*(C0 + u*(C1 + u*(C2 + u*u*C4)))
end

function lb_volume(model::KolafaNezbedaLJ,T,z)
    d̃ = d_kn(T)
    return sum(z)* d̃^3 *π/6
end

export KolafaNezbeda, KolafaNezbedaLJ