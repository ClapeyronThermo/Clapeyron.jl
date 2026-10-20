struct LJRef{T,I,LJ} <: CorrespondingStatesModel
    components::Vector{String}
    params::ChainLJScalingParam{T}
    idealmodel::I
    ljmodel::LJ
end

Base.eltype(model::LJRef) = Base.promote_eltype(model.idealmodel,model.ljmodel,model.params)


"""
    LJRef <: CorrespondingStatesModel

    LJRef(components; idealmodel = BasicIdeal, ljmodel = Thol2016LJ(), userlocations = String[], verbose = false)

## Input parameters
- `sigma`: Single Parameter (`Float64`) - particle size `[Å]`
- `epsilon`: Single Parameter (`Float64`) - dispersion energy `[K]`
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `k`: Pair Parameter (`Float64`) (optional)

- `sigma` mixing coefficient

## Model Parameters
- `sigma`: Pair Parameter (`Float64`) - particle size `[m]`
- `epsilon`: Pair Parameter (`Float64`) - dispersion energy `[K]`
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `segment`: 

## Description

vdW-1 fluid scaling applied to a generic leonard jones model:

```julia
a_res(model,V,T,z) = a_res(reference_model(model),V/h,T/f)
```

"""
LJRef

function LJRef(components;
    ljmodel = Thol2016LJ(),
    idealmodel = BasicIdeal,
    userlocations = String[],
    ideal_userlocations = String[],
    verbose = false)
    
    _components = format_components(components)
    params1 = getparams(_components, default_locations(PCSAFT); userlocations = userlocations, verbose = verbose)
    params2 = transform_params(ChainLJScaling,params1)
    packagedparams = build_eosparam(ChainLJScalingParam,params2)
    init_idealmodel = init_model(idealmodel,_components,ideal_userlocations,verbose)
    model = LJRef(_components,packagedparams,init_idealmodel, ljmodel)
    set_reference_state!(model,reference_state;verbose)
    return model
end

