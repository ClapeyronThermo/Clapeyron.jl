const vdWParam = ABCubicParam

abstract type vdWModel <: ABCubicModel end

struct vdW{T <: IdealModel,α,c,M} <: vdWModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    alpha::α
    mixing::M
    translation::c
    params::vdWParam
    idealmodel::T
    references::Array{String,1}
end

@registermodel vdW
export vdW

"""
    vdW(components::Vector{String};
    idealmodel=BasicIdeal,
    alpha = NoAlpha,
    mixing = vdW1fRule,
    activity=nothing,
    translation=NoTranslation,
    userlocations=String[], 
    ideal_userlocations=String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    verbose=false)

## Input parameters
- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Pressure `[Pa]`
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `k`: Pair Parameter (`Float64`)

## Model Parameters
- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Pressure `[Pa]`
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `a`: Pair Parameter (`Float64`)
- `b`: Pair Parameter (`Float64`)

## Input models
- `idealmodel`: Ideal Model
- `alpha`: Alpha model
- `mixing`: Mixing model
- `activity`: Activity Model, used in the creation of the mixing model.
- `translation`: Translation Model

## Description

van der Wals Equation of state.

```
P = RT/(V-Nb) + a•α(T)/V²
```

## References

1. van der Waals JD. Over de Continuiteit van den Gasen Vloeistoftoestand. PhD thesis, University of Leiden; 1873

"""
vdW

function vdW(components::Vector{String}; idealmodel=BasicIdeal,
    alpha = NoAlpha,
    mixing = vdW1fRule,
    activity=nothing,
    translation=NoTranslation,
    userlocations=String[], 
    ideal_userlocations=String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    verbose=false)
    params = getparams(components, ["properties/critical.csv", "properties/molarmass.csv","SAFT/PCSAFT/PCSAFT_unlike.csv"]; userlocations=userlocations, verbose=verbose)
    k  = params["k"]
    pc = params["pc"]
    Mw = params["Mw"]
    Tc = params["Tc"]
    init_mixing = init_model(mixing,components,activity,mixing_userlocations,activity_userlocations,verbose)
    a,b = ab_premixing(vdW,init_mixing,Tc,pc,k)
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    init_translation = init_model(translation,components,translation_userlocations,verbose)
    init_alpha = init_model(alpha,components,alpha_userlocations,verbose)
    icomponents = 1:length(components)
    packagedparams = vdWParam(a,b,Tc,pc,Mw)
    references = String[]
    model = vdW(components,icomponents,init_alpha,init_mixing,init_translation,packagedparams,init_idealmodel,references)
    return model
end

function ab_consts(::Type{<:vdWModel})
    Ωa =  27/64
    Ωb =  1/8
    return Ωa,Ωb
end

cubic_Δ(model::vdWModel,z) = (0.0,0.0)

#when either Δ1 or Δ2 is equal to zero, requires defining a_res
function a_res(model::vdWModel, V, T, z,_data = data(model,V,T,z))
    n,ā,b̄,c̄ = _data
    RT⁻¹ = 1/(R̄*T)
    ρt = (V/n+c̄)^(-1) # translated density
    ρ  = n/V
    return -log(1+(c̄-b̄)*ρ) - ā*ρt*RT⁻¹
    #
    #return -log(V-n*b̄) - ā*n/(R̄*T*V) + log(V)
end   