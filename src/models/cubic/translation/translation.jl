"""
    translation(model::CubicModel,V,T,z,translation_model::TranslationModel)  

Interface function used in cubic models. it should return a vector of cᵢ. such as `Ṽ = V + mixing(c,z)`

## Example:

```julia
function α_function(model::CubicModel,V,T,z,translation_model::RackettTranslation)
    Tc = model.params.Tc.values
    Pc = model.params.Pc.values
    Vc = translation_model.params.Vc.values
    R = Clapeyron.R̄
    Zc = Pc .* Vc ./ (R .* Tc)
    c = 0.40768 .* (0.29441 .- Zc) .* R .* Tc ./ Pc
    return c
end
```
"""
function translation end

include("NoTranslation.jl")
include("Rackett.jl")
include("Peneloux.jl")
include("MT.jl")
include("ConstantTranslation.jl")
#include("Clausius.jl")
