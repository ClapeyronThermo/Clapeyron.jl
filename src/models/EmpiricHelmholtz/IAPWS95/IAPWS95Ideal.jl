struct IAPWS95Ideal <: IdealModel
    components::Vector{String}
    references::Vector{String}
end

IAPWS95Ideal() = IAPWS95Ideal(["water"],["IAPWS R6-95(2018)"])
function IAPWS95Ideal(components;userlocations = String[],verbose = false)
    return IAPWS95Ideal()
end

idealmodel(model::IAPWS95) = IAPWS95Ideal(model.components,model.references)

"""
    IAPWS95Ideal <: IdealModel
    IAPWS95Ideal(components; 
    userlocations::Array{String,1}=String[], 
    verbose=false)

    IAPWS95Ideal()

## Input parameters

None

## Description

IAPWS95 ideal helmholtz model for use in other models. Only valid for water. Check [`IAPWS95`](@ref) for more information.

## References

1. Wagner, W., & Pruß, A. (2002). The IAPWS formulation 1995 for the thermodynamic properties of ordinary water substance for general and scientific use. Journal of physical and chemical reference data, 31(2), 387–535. [doi:10.1063/1.1461829](https://doi.org/10.1063/1.1461829)
2. IAPWS R6-95 (2018). Revised Release on the IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary Water Substance for General and Scientific Use

"""
IAPWS95Ideal

function a_ideal(model::IAPWS95Ideal,V,T,z=SA[1.0])
    Σz = only(z) #single component
    v = V/Σz
    mass_v =  v*1000.0*0.055508472036052976
    rho = one(mass_v)/mass_v
    δ = rho*0.003105590062111801 #/322
    τ = 647.096/T
    return 0.9999890238768239*iapws_f0(model,δ,τ)
end