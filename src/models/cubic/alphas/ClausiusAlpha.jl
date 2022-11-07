abstract type ClausiusAlphaModel <: AlphaModel end

struct ClausiusAlphaParam <: EoSParam
end

@newmodelsimple ClausiusAlpha ClausiusAlphaModel ClausiusAlphaParam
export ClausiusAlpha

"""
    ClausiusAlpha <: ClausiusAlphaModel
    
    ClausiusAlpha(components::Vector{String};
    userlocations::Vector{String}=String[],
    verbose::Bool=false)

## Input Parameters

- `w`: Single Parameter (`Float64`)

## Model Parameters

- `acentricfactor`: Single Parameter (`Float64`)

## Description

Cubic alpha `(α(T))` model. Default for [`Clausius`](@ref)  and [`Berthelot`]
```
αᵢ = 1/Trᵢ
Trᵢ = T/Tcᵢ
```

"""
ClausiusAlpha

function ClausiusAlpha(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false)
    packagedparams = ClausiusAlphaParam()
    model = ClausiusAlpha(packagedparams, verbose=verbose)
    return model
end

ClausiusAlpha() = ClausiusAlpha(ClausiusAlphaParam())

function α_function(model::CubicModel,V,T,z,alpha_model::ClausiusAlphaModel)
    Tc = model.params.Tc.values
    α = zeros(typeof(1.0*T),length(Tc))
    for i in @comps
        α[i] = Tc[i]/T
    end
    return α
end

function α_function(model::CubicModel,V,T,z::SingleComp,alpha_model::ClausiusAlphaModel)
    Tc = model.params.Tc.values[1]
    α = Tc/T
end

is_splittable(::ClausiusAlpha) = false