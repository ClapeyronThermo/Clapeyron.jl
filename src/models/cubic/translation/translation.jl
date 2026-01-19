"""
    translation(model::CubicModel,V,T,z,translation_model::TranslationModel)  

Interface function used in cubic models. It should return a vector of cᵢ. Such as `Ṽ = V + mixing(c,z)`

## Example:

```julia
function translation(model::CubicModel,V,T,z,translation_model::RackettTranslation)
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

translation(model::CubicModel,V,T,z) = translation(model,V,T,z,model.translation)

recombine_translation!(model::CubicModel,translation_model) = translation_model

include("NoTranslation.jl")
include("Rackett.jl")
include("Peneloux.jl")
include("MT.jl")
include("ConstantTranslation.jl")
include("VTPR.jl")
#include("Clausius.jl")
