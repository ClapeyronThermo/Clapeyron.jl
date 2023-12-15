abstract type SolidHfusModel <: EoSModel end

struct SolidHfusParam <: EoSParam
    Hfus::SingleParam{Float64}
    Tm::SingleParam{Float64}
    CpSL::SingleParam{Float64}
end

@newmodelsimple SolidHfus SolidHfusModel SolidHfusParam

"""
    SolidHfusModel <: EoSModel
    
    SolidHfus(components;
    userlocations=String[],
    verbose::Bool=false)

## Parameters

- `Hfus`: Single Parameter (`Float64`) - Enthalpy of Fusion `[J/mol]`
- `Tm`: Single Parameter (`Float64`) - Melting Temperature `[K]`
- `CpSL`: Single Parameter (`Float64`) - Heat Capacity of the Solid-Liquid Phase Transition `[J/mol/K]`

## Description

Approximation of the excess chemical potential in the solid phase (`CpSL` is not necessary by default): 
```
ln(x_iγ_i) = Hfus*T*(1/Tm-1/T)-CpSL/R̄*(Tm/T-1-log(Tm/T))
```
"""
SolidHfus
default_locations(::Type{SolidHfus}) = ["solids/fusion.csv"]
default_references(::Type{SolidHfus}) = String[]
default_ignore_missing_singleparams(::Type{SolidHfus}) = ["CpSL"]

function chemical_potential(model::SolidHfusModel, p, T, z)
    Hfus = model.params.Hfus.values
    Tm = model.params.Tm.values
    CpSL = model.params.CpSL.values
    return @. Hfus*T*(1/Tm-1/T)-CpSL/Rgas()*(Tm/T-1-log(Tm/T))
end

export SolidHfus