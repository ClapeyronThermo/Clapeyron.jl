abstract type AntoineEqSatModel <: SaturationModel end

struct AntoineEqSatParam <: EoSParam 
    Tc::SingleParam{Float64}
    Pc::SingleParam{Float64}
    A::SingleParam{Float64}
    B::SingleParam{Float64}
    C::SingleParam{Float64}
    Tmin::SingleParam{Float64}
    Tmax::SingleParam{Float64}
end

@newmodelsimple AntoineEqSat AntoineEqSatModel AntoineEqSatParam

"""
    AntoineEqSat <: SaturationModel
    
    AntoineEqSat(components;
    userlocations = String[],
    verbose::Bool=false)

## Input Parameters

- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Pressure `[Pa]`
- `A`: Single Parameter (`Float64`) - First coefficient `[dimensionless]`
- `B`: Single Parameter (`Float64`) - Second coefficent `[°C]`
- `C`: Single Parameter (`Float64`) - Third coefficent `[°C]`
- `Tmin`: Single Parameter (`Float64`)  - Mininum Temperature range `[K]`
- `Tmax`: Single Parameter (`Float64`)  - Maximum Temperature range `[K]`

## Model Parameters

- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Pressure `[Pa]`
- `A`: Single Parameter (`Float64`) - First coefficient `[dimensionless]`
- `B`: Single Parameter (`Float64`) - Second coefficent `[°C]`
- `C`: Single Parameter (`Float64`) - Third coefficent `[°C]`
- `Tmin`: Single Parameter (`Float64`)  - Mininum Temperature range `[K]`
- `Tmax`: Single Parameter (`Float64`)  - Maximum Temperature range `[K]`

## Description

Antoine Equation for saturation pressure:
```
pₛ / Pa = 10^(A + B/(T - 273.15 K + C))
```

## References

1. Antoine, C. (1888), «Tensions des vapeurs; nouvelle relation entre les tensions et les températures», Comptes Rendus des Séances de l'Académie des Sciences 107: 681-684, 778-780, 836-837.
"""
AntoineEqSat
default_locations(::Type{AntoineEqSat}) = ["properties/critical.csv","Correlations/saturation_correlations/Antoine_like.csv"]

function crit_pure(model::AntoineEqSatModel)
    single_component_check(crit_pure,model)
    tc = only(model.params.Tc.values)
    pc = only(model.params.Pc.values)
    return (tc,pc,NaN)
end

function saturation_pressure_impl(model::AntoineEqSatModel,T,::SaturationCorrelation; )
    nan = zero(T)/zero(T)
    Tc = only(model.params.Tc.values)
    A = only(model.params.A.values)
    B = only(model.params.B.values)
    C = only(model.params.C.values)
    Tmin = only(model.params.Tmin.values)
    Tmax = only(model.params.Tmax.values)

    T > Tc && (return nan,nan,nan)
    Tmin <= T <= Tmax || (return nan,nan,nan)
    psat = 10^(A - B/(T + C - 273.15))
    return psat,nan,nan
end

export AntoineEqSat
