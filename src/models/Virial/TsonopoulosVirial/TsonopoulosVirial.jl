struct TsonopoulosVirialParam <: EoSParam
    Tc::SingleParam{Float64}
    Pc::SingleParam{Float64}
    acentricfactor::SingleParam{Float64}
    Mw::SingleParam{Float64}
end

abstract type TsonopoulosVirialModel <: IdealModel end
@newmodel TsonopoulosVirial VirialModel TsonopoulosVirialParam


"""
    TsonopoulosVirial <: VirialModel
    TsonopoulosVirial(components::Array{String,1}; 
    userlocations::Array{String,1}=String[], 
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

## Description

Virial model using Corresponding State Principles
"""
TsonopoulosVirial

export TsonopoulosVirial
function TsonopoulosVirial(components;
    idealmodel=BasicIdeal,
    userlocations=String[],
    ideal_userlocations=String[],
    verbose=false)
    params = getparams(components, ["properties/molarmass"]; userlocations=userlocations, verbose=verbose)
    Mw = params["Mw"]
    Tc = params["Tc"]
    Pc = params["pc"]
    acentricfactor = params["w"]
    packagedparams = TsonopoulosVirialParam(Tc,Pc,acentricfactor,Mw)
    references = String[]
    return PCSAFT(packagedparams, idealmodel; ideal_userlocations, references, verbose)
end

function second_virial_coefficient(model::TsonopoulosVirialModel,T,z=SA[1.0])
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
            Tcj = Tc[i]
            Pcj = Pc[i]
            ωj = ω[j]
            #mixing, todo: implement mixing at the EoS level, like cubics
            Tcij = sqrt(Tci*Tcj)
            Pcij = 0.5*(Pci+Pcj)
            ωij =  0.5*(ωi+ωj)
            zi = z[i]
            τ = Tci/T
            B0 = 0.335 - 0.43*τ -  0.1385*τ*τ - 0.0121*τ*τ*τ - 0.000607*τ^8
            B1 = 0.0637 + 0.331*τ*τ - 0.423*τ*τ*τ - 0.008*τ^8
            B +=  2*zi*zj*(B0+ ωij*B1)*R̄*Tcij/Pcij
        end
    end
    return B/sum(z)
end