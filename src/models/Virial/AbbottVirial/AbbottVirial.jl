struct AbbottVirialParam <: EoSParam
    Tc::SingleParam{Float64}
    Pc::SingleParam{Float64}
    acentricfactor::SingleParam{Float64}
    Mw::SingleParam{Float64}
end

@newmodel AbbottVirial SecondVirialModel AbbottVirialParam false
default_locations(::Type{AbbottVirial}) = ["properties/critical.csv", "properties/molarmass.csv"]

"""
    AbbottVirial <: SecondVirialModel
    AbbottVirial(components;
                idealmodel = BasicIdeal,
                userlocations = String[],
                ideal_userlocations = String[],
                verbose = false)

## Input parameters

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
B₀ = 0.083 + 0.422/Trᵢⱼ^1.6
B₁ = 0.139 - 0.172/Trᵢⱼ^4.2
Trᵢⱼ = T/Tcᵢⱼ
Tcᵢⱼ = √TcᵢTcⱼ
Pcᵢⱼ = (Pcᵢ + Pcⱼ)/2
ωᵢⱼ = (ωᵢ + ωⱼ)/2
```

## Model Construction Examples
```julia
# Using the default database
model = AbbottVirial("water") #single input
model = AbbottVirial(["water","ethanol"]) #multiple components
model = AbbottVirial(["water","ethanol"], idealmodel = ReidIdeal) #modifying ideal model

# Passing a prebuilt model

my_idealmodel = MonomerIdeal(["neon","hydrogen"];userlocations = (;Mw = [20.17, 2.]))
model = AbbottVirial(["neon","hydrogen"],idealmodel = my_idealmodel)

# User-provided parameters, passing files or folders
model = AbbottVirial(["neon","hydrogen"]; userlocations = ["path/to/my/db","critical.csv"])

# User-provided parameters, passing parameters directly

model = AbbottVirial(["neon","hydrogen"];
        userlocations = (;Tc = [44.492,33.19],
                        Pc = [2679000, 1296400],
                        Mw = [20.17, 2.],
                        acentricfactor = [-0.03,-0.21])
                    )
```

## References
1. Smith, H. C. Van Ness Joseph M. Introduction to Chemical Engineering Thermodynamics 4E 1987.
"""
AbbottVirial

export AbbottVirial

function second_virial_coefficient_impl(model::AbbottVirial,T,z=SA[1.0])
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
        B0 = 0.083 - 0.422*τ^1.6
        B1 = 0.139 + 0.179*τ^4.2
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
            B0 = 0.083 - 0.422*τ^1.6
            B1 = 0.139 + 0.179*τ^4.2
            B +=  2*zi*zj*(B0+ ωij*B1)*R̄*Tcij/Pcij
        end
    end
    return B/sum(z)
end