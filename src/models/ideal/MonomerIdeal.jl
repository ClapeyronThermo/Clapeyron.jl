struct MonomerIdealParam <: EoSParam
    Mw::SingleParam{Float64}
end

abstract type MonomerIdealModel <: IdealModel end
@newmodelsimple MonomerIdeal MonomerIdealModel MonomerIdealParam


"""
    MonomerIdeal <: MonomerIdealModel
    MonomerIdeal(components::Array{String,1}; 
    userlocations::Array{String,1}=String[], 
    verbose=false)

## Input parameters

- `Mw`: Single Parameter (`Float64`)

## Model Parameters

None

Monomer Ideal Model, calculated as a function of the thermal wavelength `Λ`
```
    Λᵢ = h/√(kᵦTMwᵢ/Nₐ)
    a₀ = A₀/nRT = ∑xᵢlog((Λᵢ^3)*n*xᵢ*Nₐ/V)
```
"""

export MonomerIdeal
function MonomerIdeal(components::Array{String,1}; userlocations::Array{String,1}=String[], verbose=false)
    params = getparams(components, ["properties/molarmass"]; userlocations=userlocations, verbose=verbose)
    Mw = params["Mw"]
    packagedparams = MonomerIdealParam(Mw)
    return MonomerIdeal(packagedparams)
end

function a_ideal(model::MonomerIdealModel, v, T, z)
    Mw = model.params.Mw.values
    res = zero(V+T+first(z))
    for i in @comps
        Λᵢ = h/√(k_B*T*Mw[i]/N_A)
        res +=  z[i]*log(z[i]*N_A/v*Λᵢ^3)
    end
    return res/sum(z) - 1
end
