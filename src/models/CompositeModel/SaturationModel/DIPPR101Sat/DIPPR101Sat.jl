abstract type DIPPR101SatModel <: SaturationModel end

struct DIPPR101SatParam <: EoSParam 
    Tc::SingleParam{Float64}
    Pc::SingleParam{Float64}
    A::SingleParam{Float64}
    B::SingleParam{Float64}
    C::SingleParam{Float64}
    D::SingleParam{Float64}
    E::SingleParam{Float64}
    Tmin::SingleParam{Float64}
    Tmax::SingleParam{Float64}
end

@newmodelsimple DIPPR101Sat DIPPR101SatModel DIPPR101SatParam

"""
    DIPPR101Sat <: SaturationModel
    
    DIPPR101Sat(components;
    userlocations = String[],
    verbose::Bool=false)

## Input Parameters

- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Pressure `[Pa]`
- `A`: Single Parameter (`Float64`)
- `B`: Single Parameter (`Float64`)
- `C`: Single Parameter (`Float64`)
- `D`: Single Parameter (`Float64`)
- `E`: Single Parameter (`Float64`)
- `Tmin`: Single Parameter (`Float64`)  - mininum Temperature range `[K]`
- `Tmax`: Single Parameter (`Float64`)  - maximum Temperature range `[K]`

## Model Parameters

- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Pressure `[Pa]`
- `A`: Single Parameter (`Float64`)
- `B`: Single Parameter (`Float64`)
- `C`: Single Parameter (`Float64`)
- `D`: Single Parameter (`Float64`)
- `E`: Single Parameter (`Float64`)
- `Tmin`: Single Parameter (`Float64`)  - mininum Temperature range `[K]`
- `Tmax`: Single Parameter (`Float64`)  - maximum Temperature range `[K]`

## Description

DIPPR 101 Equation for saturation pressure:
```
psat(T) = exp(A + B/T + C•log(T) + D•T^E)
```

## References

1. Design Institute for Physical Properties, 1996. DIPPR Project 801 DIPPR/AIChE
"""
DIPPR101Sat
default_locations(::Type{DIPPR101Sat}) = ["properties/critical.csv","Correlations/saturation_correlations/dippr101_like.csv"]

function crit_pure(model::DIPPR101SatModel)
    single_component_check(crit_pure,model)
    tc = only(model.params.Tc.values)
    pc = only(model.params.Pc.values)
    return (tc,pc,NaN)
end

function saturation_pressure_impl(model::DIPPR101SatModel,T,method::SaturationCorrelation)
    nan = zero(T)/zero(T)
    tc = only(model.params.Tc.values)
    A = only(model.params.A.values)
    B = only(model.params.B.values)
    C = only(model.params.C.values)
    D = only(model.params.D.values)
    E = only(model.params.E.values)
    Tmin = only(model.params.Tmin.values)
    Tmax = only(model.params.Tmax.values)

    T > tc && (return nan,nan,nan)
    Tmin <= T <= Tmax || (return nan,nan,nan)
    psat = exp(A + B/T + C*log(T) + D*T^E)
    return psat,nan,nan
end

export DIPPR101Sat
