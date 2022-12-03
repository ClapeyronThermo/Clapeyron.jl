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

- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`

## Model Parameters

None

## Description

Monomer Ideal Model, result obtained from statistical mechanics `Λ`
```
    Λᵢ = h/√(kᵦTMwᵢ/Nₐ)    
    a₀ = A₀/nRT = ∑xᵢlog(ρᵢΛᵢ^3)
```
"""
MonomerIdeal

export MonomerIdeal
function MonomerIdeal(components::Array{String,1}; userlocations::Array{String,1}=String[], verbose=false)
    params = getparams(components, ["properties/molarmass.csv"]; userlocations=userlocations, verbose=verbose)
    Mw = params["Mw"]
    packagedparams = MonomerIdealParam(Mw)
    return MonomerIdeal(packagedparams)
end

recombine_impl!(model::MonomerIdealModel) = model


function a_ideal(model::MonomerIdealModel, V, T, z)
    Mw = model.params.Mw.values
    res = zero(V+T+first(z))
    for i in @comps
        Mwᵢ = Mw[i]*0.001
        Λᵢ = h/√(k_B*T*Mwᵢ/N_A)
        res +=  z[i]*log(z[i]*N_A/V*Λᵢ^3)
    end
    return res/sum(z) - 1
end

molecular_weight(model::MonomerIdealModel,z) = comp_molecular_weight(mw(model),z)