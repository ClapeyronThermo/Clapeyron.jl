abstract type MathiasCopemanAlphaModel <: AlphaModel end

struct MathiasCopemanAlphaParam{T} <: ParametricEoSParam{T}
    c1::SingleParam{T}
    c2::SingleParam{T}
    c3::SingleParam{T}
end

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
2.  S. Horstmann, A. Jabłoniec, J. Krafczyk, K. Fischer, and J. Gmehling: PSRK Group Contribution Equation of State: Comprehensive Revision and Extension IV, Including Critical Constants and α-Function Parameters for 1000 Components, Fluid Phase Equilibria 227 (2005) 157–164, DOI: https://doi.org/10.1016/j.fluid.2004.11.002.
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

function α_function(model::CubicModel,V,T,z,alpha_model::MathiasCopemanAlpha)
    Tc = model.params.Tc.values
    c1 = alpha_model.params.c1.values
    c2 = alpha_model.params.c2.values
    c3 = alpha_model.params.c3.values
    α = zeros(typeof(T*1.0),length(Tc))
    for i in @comps
        Tr = T/Tc[i]
        _1_Tr = 1 - √(Tr)
        α[i] = (1 + c1[i]*_1_Tr + c2[i]*_1_Tr^2 + c3[i]*_1_Tr^3)^2
    end
    return α
end
