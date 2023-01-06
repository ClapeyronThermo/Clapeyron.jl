struct TsonopoulosVirialParam <: EoSParam
    Tc::SingleParam{Float64}
    Pc::SingleParam{Float64}
    acentricfactor::SingleParam{Float64}
    Mw::SingleParam{Float64}
end

@newmodel TsonopoulosVirial SecondVirialModel TsonopoulosVirialParam

"""
    TsonopoulosVirial <: SecondVirialModel
    TsonopoulosVirial(components;
            idealmodel=BasicIdeal,
            userlocations=String[],
            ideal_userlocations=String[],
            verbose=false)

## Input parameters

- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Pressure `[Pa]`
- `w`: Single Parameter (`Float64`) - Acentric Factor
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`

## Model Parameters

- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Pressure `[Pa]`
- `acentricfactor`: Single Parameter (`Float64`) - Acentric Factor
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`

## Input models
- `idealmodel`: Ideal Model

## Description

Virial model using Corresponding State Principles:
```
B = ∑xᵢxⱼBᵢⱼ
Bᵢⱼ = BrᵢⱼRTcᵢⱼ/Pcᵢⱼ
Brᵢⱼ = B₀ + ωᵢⱼB₁
B₀ = 0.1445 - 0.330/Trᵢⱼ - 0.1385/Trᵢⱼ^2 - 0.0121/Trᵢⱼ^3 - 0.000607/Trᵢⱼ^8
B₁ = 0.0637 + 0.331/Trᵢⱼ - 0.423/Trᵢⱼ^2 - 0.423/Trᵢⱼ^3 - 0.008/Trᵢⱼ^8
Trᵢⱼ = T/Tcᵢⱼ
Tcᵢⱼ = √TcᵢTcⱼ
Pcᵢⱼ = (Pcᵢ + Pcⱼ)/2
ωᵢⱼ = (ωᵢ + ωⱼ)/2
```

## References

1. Tsonopoulos, C. (1974). An empirical correlation of second virial coefficients. AIChE Journal. American Institute of Chemical Engineers, 20(2), 263–272. [doi:10.1002/aic.690200209](https://doi.org/10.1002/aic.690200209)

"""
TsonopoulosVirial

export TsonopoulosVirial
function TsonopoulosVirial(components;
    idealmodel=BasicIdeal,
    userlocations=String[],
    ideal_userlocations=String[],
    verbose=false)
    params = getparams(components,["properties/critical.csv", "properties/molarmass.csv"]; userlocations=userlocations, verbose=verbose)
    Mw = params["Mw"]
    Tc = params["Tc"]
    Pc = params["Pc"]
    acentricfactor = params["acentricfactor"]
    packagedparams = TsonopoulosVirialParam(Tc,Pc,acentricfactor,Mw)
    references = String["10.1002/aic.690200209"]
    return TsonopoulosVirial(packagedparams, idealmodel; ideal_userlocations, references, verbose)
end

function second_virial_coefficient_impl(model::TsonopoulosVirial,T,z=SA[1.0])
    B = zero(T+first(z))
    Tc = model.params.Tc.values
    ω = model.params.acentricfactor.values
    Pc = model.params.Pc.values
    for i in 1:length(z)
        Tci = Tc[i]
        Pci = Pc[i]
        ωi = ω[i]
        zi = z[i]
        τ = Tci/T
        B0 = 0.335 - 0.43*τ -  0.1385*τ*τ - 0.0121*τ*τ*τ - 0.000607*τ^8
        B1 = 0.0637 + 0.331*τ*τ - 0.423*τ*τ*τ - 0.008*τ^8
        B +=  zi*zi*(B0+ ωi*B1)*R̄*Tci/Pci
        for j in 1:i-1
            Tcj = Tc[j]
            Pcj = Pc[j]
            ωj = ω[j]
            #mixing, todo: implement mixing at the EoS level, like cubics
            Tcij = sqrt(Tci*Tcj)
            Pcij = 0.5*(Pci+Pcj)
            ωij =  0.5*(ωi+ωj)
            zj = z[j]
            τ = Tci/T
            B0 = 0.335 - 0.43*τ -  0.1385*τ*τ - 0.0121*τ*τ*τ - 0.000607*τ^8
            B1 = 0.0637 + 0.331*τ*τ - 0.423*τ*τ*τ - 0.008*τ^8
            B +=  2*zi*zj*(B0+ ωij*B1)*R̄*Tcij/Pcij
        end
    end
    return B/sum(z)
end