abstract type LeeKeslerSatModel <: SaturationModel end

struct LeeKeslerSatParam <: EoSParam 
    Tc::SingleParam{Float64}
    Pc::SingleParam{Float64}
    acentricfactor::SingleParam{Float64}
end

@newmodelsimple LeeKeslerSat LeeKeslerSatModel LeeKeslerSatParam

"""
    LeeKeslerSat <: SaturationModel
    
    LeeKeslerSat(components::Vector{String};
    userlocations::Vector{String}=String[],
    verbose::Bool=false)

## Input Parameters

- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Pressure `[Pa]`
- `acentricfactor`: Single Parameter (`Float64`) - acentric factor

## Description

Lee-Kesler correlation for saturation pressure:
```
psat(T) = f₀ + ω•f₁
Tr = T/Tc
f₀ = 5.92714 - 6.09648/Tr - 1.28862•log(Tr) + 0.169347•Tr⁶
f₁ = 15.2518 - 15.6875/Tr - 13.4721•log(Tr) + 0.43577•Tr⁶
```

## References

1. Lee, B. I., & Kesler, M. G. (1975). A generalized thermodynamic correlation based on three-parameter corresponding states. AIChE journal. American Institute of Chemical Engineers, 21(3), 510–527. [doi:10.1002/aic.690210313](https://doi.org/10.1002/aic.690210313)
"""
LeeKeslerSat

function LeeKeslerSat(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false)
    params = getparams(components, ["properties/critical.csv"]; userlocations=userlocations, verbose=verbose)
    acentricfactor = params["acentricfactor"]
    Tc = params["Tc"]
    Pc = params["Pc"]
    packagedparams = LeeKeslerSatParam(Tc,Pc,acentricfactor)
    model = LeeKeslerSat(packagedparams, verbose=verbose)
    return model
end 

function crit_pure(model::LeeKeslerSatModel)
    tc = only(model.params.Tc.values)
    pc = only(model.params.Pc.values)
    return (tc,pc,NaN)
end

function saturation_pressure_impl(model::LeeKeslerSatModel,T,method::SaturationCorrelation)
    nan = zero(T)/zero(T)
    ω = only(model.params.acentricfactor.values)
    tc = only(model.params.Tc.values)
    pc = only(model.params.Pc.values)
    T > tc && (return nan,nan,nan)
    
    tr = T/tc
    trinv = inv(tr)
    lntr = log(tr)
    tr6 = tr^6
    
    f0 = 5.92714 - 6.09648*trinv - 1.28862*lntr + 0.169347*tr6
    f1 = 15.2518 - 15.6875*trinv - 13.4721*lntr + 0.43577*tr6
    lnpr = f0 + ω*f1
    psat = exp(lnpr)*pc
    return psat,nan,nan
end

export LeeKeslerSat
