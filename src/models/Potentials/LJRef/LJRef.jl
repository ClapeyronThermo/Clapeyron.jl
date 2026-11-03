struct LJRef{T,I,LJ} <: CorrespondingStatesModel
    components::Vector{String}
    params::ChainLJScalingParam{T}
    idealmodel::I
    ljmodel::LJ
    molecular_units::Bool
end

Base.eltype(model::LJRef) = Base.promote_eltype(model.idealmodel,model.ljmodel,model.params)


"""
    LJRef <: CorrespondingStatesModel

    LJRef(components;
        idealmodel = BasicIdeal,
        ljmodel = Thol2016LJ(),
        userlocations = String[],
        molecular_units = false,
        verbose = false)

## Input parameters
- `sigma`: Single Parameter (`Float64`) - particle size `[Å]`
- `epsilon`: Single Parameter (`Float64`) - dispersion energy `[K]`
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `k`: Pair Parameter (`Float64`) (optional)

## Model Parameters
- `sigma`: Pair Parameter (`Float64`) - particle size `[m]`
- `epsilon`: Pair Parameter (`Float64`) - dispersion energy `[K]`
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`

## Description

vdW-1 fluid scaling applied to a generic leonard jones model.

```julia
a_res(model,V,T,z) = a_res(reference_model(model),V/h,T/f)
```

The `molecular_units` argument sets if the scaling is done in SI units or molecular units.
On SI units, the input sigma diameter is scaled from `[Å]` to `[m]` (once, at construction time), the volume is scaled by the avogadro number, and the gas constant is set to `R = 8.314...`;
on molecular units, no additional scaling is applied and `R = 1.0`
"""
function LJRef(components;
    ljmodel::EoSModel = Thol2016LJ(),
    idealmodel = BasicIdeal,
    userlocations = String[],
    ideal_userlocations = String[],
    molecular_units = false,
    reference_state = nothing,
    verbose = false)

    _components = format_components(components)
    params1 = getparams(_components, default_locations(PCSAFT); userlocations = userlocations, verbose = verbose)
    params2 = transform_params_chainljscaling(params1,_components,molecular_units)
    packagedparams = build_eosparam(ChainLJScalingParam,params2)
    T = eltype(packagedparams)
    init_idealmodel = init_model(idealmodel,_components,ideal_userlocations,verbose)
    model = LJRef{T,typeof(init_idealmodel),typeof(ljmodel)}(_components, packagedparams, init_idealmodel, ljmodel, molecular_units)
    set_reference_state!(model,reference_state;verbose)
    return model
end

reference_model(model::LJRef) = model.ljmodel
scaling_model(model::LJRef) = ChainLJScaling(model.params,model.molecular_units)

function recombine_impl!(model::LJRef)
    recombine!(model.ljmodel)
    recombine_lj_cs!(model,nothing,nothing)
end

set_k!(model::LJRef,k) = recombine_lj_cs!(model,k,nothing)
set_l!(model::LJRef,l) = recombine_lj_cs!(model,nothing,l)

export LJRef