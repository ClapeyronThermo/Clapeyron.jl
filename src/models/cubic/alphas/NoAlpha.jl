abstract type NoAlphaModel <: AlphaModel end

struct NoAlphaParam <: EoSParam
end

@newmodelsimple NoAlpha NoAlphaModel NoAlphaParam
export NoAlpha

"""
    NoAlpha(args...) <: NoAlphaModel
## Input Parameters
None
## Description
Cubic alpha `(α(T))` model. Default for [`vdW`](@ref) EoS
```
αᵢ = 1 ∀ i
```
"""
NoAlpha

function α_function(model::CubicModel,V,T,z,alpha_model::NoAlphaModel)
   return FillArrays.Ones{Float64}(length(z))
end