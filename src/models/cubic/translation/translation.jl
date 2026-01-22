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

function translation(model,V,T,z,translation)
    ∑z = sum(z)
    v = V/∑z
    _1 = oneunit(∑z)
    n = length(model)
    z1 = FillArrays.OneElement(_1,1,n)
    c1 = translation2(model,v,T,z1,translation,nothing,nothing,nothing)
    c = similar(z,typeof(c1))
    c[1] = c1
    for i in 1:n
        zi = FillArrays.OneElement(_1,i,n)
        ci = translation2(model,v,T,zi,translation,nothing,nothing,nothing)
        c[i] = ci
    end
    return c
end

recombine_translation!(model::CubicModel,translation_model) = translation_model

include("NoTranslation.jl")
include("Peneloux.jl")
include("MT.jl")
include("ConstantTranslation.jl")
#include("Clausius.jl")
