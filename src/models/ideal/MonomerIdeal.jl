struct MonomerIdealParam <: EoSParam
    Mw::SingleParam{Float64}
    reference_state::ReferenceState
end

abstract type MonomerIdealModel <: IdealModel end
@newmodelsimple MonomerIdeal MonomerIdealModel MonomerIdealParam

"""
    MonomerIdeal <: MonomerIdealModel

    MonomerIdeal(components; 
    userlocations = String[],
    reference_state = nothing,
    verbose = false)

## Input parameters

- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`

## Model Parameters

None

## Description

Monomer Ideal Model, result obtained from statistical mechanics `Λ`
```
    Λᵢ = h/√(kᵦTMwᵢ/Nₐ)    
    a₀ = A₀/nRT = ∑xᵢlog(ρᵢΛᵢ^3)
```

## Model Construction Examples
```
# Using the default database
idealmodel = MonomerIdeal("water") #single input
idealmodel = MonomerIdeal(["water","ethanol"]) #multiple components

# Using user-provided parameters

# Passing files or folders
idealmodel = MonomerIdeal(["neon","hydrogen"]; userlocations = ["path/to/my/db","mw.csv"])

# Passing parameters directly
idealmodel = MonomerIdeal(["neon","hydrogen"];userlocations = (;Mw = [20.17, 2.]))
```
"""
MonomerIdeal

export MonomerIdeal
default_locations(::Type{MonomerIdeal}) = mw_data()

recombine_impl!(model::MonomerIdealModel) = model

function a_ideal(model::MonomerIdealModel, V, T, z)
    Mw = model.params.Mw.values
    res = zero(V+T+first(z))
    for i in @comps
        Mwᵢ = Mw[i]*0.001
        Λᵢ = h/√(k_B*T*Mwᵢ/N_A)
        res += xlogx(z[i],N_A/V*Λᵢ^3)
    end
    return res/sum(z) - 1
end
