
struct LinearMixing <: MixingRule end

is_splittable(::LinearMixing) = false

"""
    QuadraticDeparture <: MultiFluidDepartureModel
    QuadraticDeparture(components; 
    userlocations=String[],
    verbose=false)

## Input parameters
- `beta_v`: Pair Parameter (`Float64`) - binary interaction parameter  (no units)
- `gamma_v`: Pair Parameter (`Float64`) - binary interaction parameter  (no units)
- `beta_T`: Pair Parameter (`Float64`) - binary interaction parameter  (no units)
- `gamma_T`: Pair Parameter (`Float64`) - binary interaction parameter  (no units)

## Description
Linear mixing rule for MultiParameter EoS models:

```
τ = T̄/T
δ = V̄/V
V̄ = ∑xᵢVcⱼ
T̄ = ∑xᵢTcᵢ
```
"""
function LinearMixing(components;userlocations = String[],verbose = false)
    LinearMixing()
end

function v_scale(model::MultiFluid,z,mixing::LinearMixing,∑z)
    Vc = model.params.Vc.values
    dot(Vc,z)/∑z
end

function T_scale(model::MultiFluid,z,mixing::LinearMixing,∑z)
    Tc = model.params.Tc.values  
    return dot(Tc,z)/∑z
end

export LinearMixing