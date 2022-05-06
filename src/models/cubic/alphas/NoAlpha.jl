abstract type NoAlphaModel <: AlphaModel end

BasicIdeal_SETUP = ModelOptions(
        :NoAlpha;
        supertype=NoAlphaModel,
        has_components=false,
        has_params=false,
    )

createmodel(BasicIdeal_SETUP; verbose=true)

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

is_splittable(::NoAlpha) = false
