abstract type MathiasCopemanAlphaModel <: AlphaModel end

struct MathiasCopemanAlphaParam{T} <: ParametricEoSParam{T}
    c1::SingleParam{T}
    c2::SingleParam{T}
    c3::SingleParam{T}
end

MathiasCopemanAlphaParam(c1,c2,c3) = build_parametric_param(MathiasCopemanAlphaParam,c1,c2,c3)

@newmodelsimple MathiasCopemanAlpha MathiasCopemanAlphaModel MathiasCopemanAlphaParam
export MathiasCopemanAlpha

"""
    MathiasCopemanAlpha <: MathiasCopemanAlphaModel
    
    MathiasCopemanAlpha(components;
    userlocations = String[],
    verbose::Bool=false)

## Input Parameters

- `c1`: Single Parameter
- `c2`: Single Parameter
- `c3`: Single Parameter

## Description

Cubic alpha `(α(T))` model using Mathias-Copeman parameters.
```
αᵢ = [1 + c₁ᵢ(1-√(Trᵢ)) + c₃ᵢ(1-√(Trᵢ))² + c₃ᵢ(1-√(Trᵢ))³)²
Trᵢ = T/Tcᵢ
```

## Model Construction Examples
```
# Passing files or folders
alpha = MathiasCopemanAlpha(["neon","hydrogen"]; userlocations = ["path/to/my/db",])

# Passing parameters directly
alpha = MathiasCopemanAlpha(["neon","hydrogen"]; userlocations=(;c1=[-0.03,-0.21], c2=[0.,0.], c3=[0.,0.]))
```

## References

1.  P. M. Mathias and T. W. Copeman: Extension of the Peng-Robinson Equation of State to Complex Mixtures: Evaluation of the Various Forms of the Local Composition Concept, Fluid Phase Equilibria 13 (1983) 91–108, DOI: https://doi.org/10.1016/0378-3812(83)80084-3.
2.  S. Horstmann, A. Jabłoniec, J. Krafczyk, K. Fischer, and J. Gmehling: PSRK Group Contribution Equation of State: Comprehensive Revision and Extension IV, Including Critical Constants and α-Function Parameters for 1000 Components, Fluid Phase Equilibria 227 (2005) 157–164, [doi:10.1016/j.fluid.2004.11.002](https://doi.org/10.1016/j.fluid.2004.11.002)
"""
MathiasCopemanAlpha

default_locations(::Type{MathiasCopemanAlpha}) = String[]
function default_ignore_missing_singleparams(::Type{T}) where T <: MathiasCopemanAlpha
    if hasfield(T,:params)
        P = fieldtype(T,:params)
        if hasfield(P,:Vc)
            return String[]
        else
            return ["Vc"]
        end
    end
    return String[]
end

function α_function(model::CubicModel,V,T,z,alpha_model::MathiasCopemanAlphaModel)
    Tc = model.params.Tc.values
    c1 = alpha_model.params.c1.values
    c2 = alpha_model.params.c2.values
    c3 = alpha_model.params.c3.values
    _1 = oneunit(eltype(alpha_model))
    α = zeros(Base.promote_eltype(Tc,T,1.0),length(Tc))
    for i in @comps
        Tr = T/Tc[i]
        _1_Tr = 1 - sqrt(Tr)
        α[i] = evalpoly(_1_Tr,(_1,c1[i],c2[i],c3[i]))^2
    end
    return α
end

function α_function(model::CubicModel,V,T,z::SingleComp,alpha_model::MathiasCopemanAlphaModel)
    Tc = model.params.Tc.values
    c1 = alpha_model.params.c1.values
    c2 = alpha_model.params.c2.values
    c3 = alpha_model.params.c3.values
    _1 = oneunit(eltype(alpha_model))
    _1_Tr = 1 - sqrt(T/Tc[1])
    return evalpoly(_1_Tr,(_1,c1[1],c2[1],c3[1]))^2
end