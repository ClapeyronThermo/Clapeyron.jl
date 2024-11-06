
struct ConstantCpLiquidParam <: EoSParam
    Cp::SingleParam{Float64} 
    reference_state::ReferenceState
end

@newmodelsimple ConstantCpLiquid PolynomialCpLiquidModel ConstantCpLiquidParam

"""
    ConstantCpLiquid <: LiquidCpModel

    ConstantCpLiquid(components; 
    userlocations = String[],
    reference_state = nothing,
    verbose = false)

## Input parameters

- `Cp`: Single Parameter (`Float64`) - isobaric heat capacity of fluid

## Model parameters

- `coeffs`: Single Parameter (`NTuple{5,Float64}`)

## Description

Incompressible fluid caloric properties. Helmholtz energy obtained via integration of heat capacity:

```
Cp = ∑Cpᵢ
```

## Model Construction Examples
```
# Using the default database
model = ConstantCpLiquid("water") #single input
model = ConstantCpLiquid(["water","ethanol"]) #multiple components

# Using user-provided parameters

# Passing files or folders
idealmodel = ConstantCpLiquid(["neon","hydrogen"]; userlocations = ["path/to/my/db","cpl.csv"])

# Passing parameters directly
idealmodel = ConstantCpLiquid(["water","butane"];
            userlocations = (Cp = [75.34, 9.487])
```
"""
ConstantCpLiquid
caloric_coefficients(model::ConstantCpLiquid) = model.params.Cp.values
export ConstantCpLiquid

