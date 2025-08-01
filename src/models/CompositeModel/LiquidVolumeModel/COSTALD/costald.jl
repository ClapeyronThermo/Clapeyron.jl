abstract type COSTALDModel <: LiquidVolumeModel end

struct COSTALDParam <: EoSParam
    Tc::SingleParam{Float64}
    Vc::SingleParam{Float64}
    acentricfactor::SingleParam{Float64}
end

@newmodelsimple COSTALD COSTALDModel COSTALDParam
"""
    COSTALD(components; 
                userlocations::Vector{String}=String[], 
                verbose::Bool=false)

## Input parameters

- `Tc`: Single Parameter (Float64) - Critical Temperature `[K]`
- `Vc`: Single Parameter (`Float64`) - Critical Volume `[m³·mol⁻¹]`
- `acentricfactor`: Single Parameter (`Float64`) - Acentric Factor

## Description

COSTALD Equation of State for saturated liquids. It is independent of the pressure.

## Model Construction Examples
```julia
# Using the default database
model = COSTALD("water") #single input
model = COSTALD(["water","ethanol"]) #multiple components

# User-provided parameters, passing files or folders
model = COSTALD(["neon","hydrogen"]; userlocations = ["path/to/my/db","critical.csv"])

# User-provided parameters, passing parameters directly

model = COSTALD(["neon","hydrogen"];
        userlocations = (;Tc = [44.492,33.19],
                        Vc = [4.25e-5, 6.43e-5],
                        acentricfactor = [-0.03,-0.21])
                    )
```

## References
Hankinson, R. W., & Thomson, G. H. (1979). A new correlation for saturated densities of liquids and their mixtures. AIChE Journal. American Institute of Chemical Engineers, 25(4), 653–663. [doi:10.1002/aic.690250412](https://doi.org/doi:10.1002/aic.690250412)
"""
COSTALD
default_locations(::Type{COSTALD}) = critical_data()
default_references(::Type{COSTALD}) = ["10.1002/aic.690250412"]

function volume_impl(model::COSTALDModel,p,T,z,phase,threaded,vol0)
    Tci = model.params.Tc.values
    Vci = model.params.Vc.values
    ωi  = model.params.acentricfactor.values
        
    Vc1 = zero(eltype(z))
    Vc23 = zero(eltype(z))
    Vc13 = zero(eltype(z))
    ω = zero(eltype(z))  
    Tc = zero(eltype(z))
    checkbounds(Tci,length(z))
    for i ∈ @comps
        zi = z[i]
        ω += zi*ωi[i]
        Vcii = Vci[i]
        Vcii13 = cbrt(Vcii)
        Vcii23 = Vcii13*Vcii13
        Vc1 += zi*Vcii
        Vc13 += zi*Vcii13
        Vc23 += zi*Vcii23
        VTi = Vcii*Tci[i]
        Tc += zi*zi*VTi
        for j ∈ 1:i-1
            VTj = Vci[j]*Tci[j]
            Tc += 2*zi*z[j]*sqrt(VTi*VTj)
        end
    end

    Vc = 0.25*(Vc1 + 3*Vc13*Vc23)
    Tc = Tc/Vc

    Tr = T/Tc
    τ = 1.0 - Tr
    τcbrt = cbrt(τ)
    Vδ = evalpoly(Tr,(-0.296123,0.386914,-0.0427258,-0.0480645))/(Tr - 1.00001)
    V0 = evalpoly(τcbrt,(1.0, -1.52816,1.43907,-0.81446,0.190454))
    return Vc*V0*(1.0 - ω*Vδ)
end

function volume_impl(model::COSTALDModel,p,T,z::SingleComp,phase,threaded,vol0)
    Tc = model.params.Tc.values |> only
    Vc = model.params.Vc.values |> only
    ω  = model.params.acentricfactor.values |> only

    Tr = T/Tc
    τ = 1.0 - Tr
    τcbrt = cbrt(τ)
    Vδ = evalpoly(Tr,(-0.296123,0.386914,-0.0427258,-0.0480645))/(Tr - 1.00001)
    V0 = evalpoly(τcbrt,(1.0, -1.52816,1.43907,-0.81446,0.190454))
    return Vc*V0*(1.0 - ω*Vδ)
end

export COSTALD