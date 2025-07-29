abstract type TwuAlphaModel <: AlphaModel end

struct TwuAlphaParam <: EoSParam
    M::SingleParam{Float64}
    N::SingleParam{Float64}
    L::SingleParam{Float64}
end

@newmodelsimple TwuAlpha TwuAlphaModel TwuAlphaParam
default_locations(::Type{TwuAlpha}) = ["alpha/Twu/Twu_like.csv"]
default_references(::Type{TwuAlpha}) = ["10.1016/0378-3812(80)80003-3"]
export TwuAlpha

"""
    TwuAlpha <: TwuAlphaModel
    Twu91Alpha = TwuAlpha
    TwuAlpha(components;
    userlocations = String[],
    verbose::Bool=false)

## Input Parameters
- `M`: Single Parameter
- `N`: Single Parameter
- `L`: Single Parameter

## Description
Cubic alpha `(α(T))` model. Default for [`VTPR`](@ref) EoS. Also known as Twu-91 alpha
```
αᵢ = Trᵢ^(N*(M-1))*exp(L*(1-Trᵢ^(N*M))
Trᵢ = T/Tcᵢ
```

## Model Construction Examples
```
# Using the default database
alpha = TwuAlpha("water") #single input
alpha = Twu91Alpha("water") #same function
alpha = TwuAlpha(["water","ethanol"]) #multiple components

# Using user-provided parameters

# Passing files or folders
alpha = TwuAlpha(["neon","hydrogen"]; userlocations = ["path/to/my/db","twu.csv"])

# Passing parameters directly
alpha = TwuAlpha(["neon","hydrogen"];
    userlocations = (;L = [0.40453, 156.21],
                    M = [0.95861, -0.0062072],
                    N = [0.8396, 5.047])
                )
```

## References
1. Twu, C. H., Lee, L. L., & Starling, K. E. (1980). Improved analytical representation of argon thermodynamic behavior. Fluid Phase Equilibria, 4(1–2), 35–44. [doi:10.1016/0378-3812(80)80003-3](https://doi.org/10.1016/0378-3812(80)80003-3)
"""
TwuAlpha

const Twu91Alpha = TwuAlpha

"""
    Twu88Alpha::TwuAlpha

    Twu88Alpha(components::Vector{String};
    userlocations = String[],
    verbose::Bool=false)

## Input Parameters
- `M`: Single Parameter
- `N`: Single Parameter (optional)
- `L`: Single Parameter

## Model Parameters
- `M`: Single Parameter
- `N`: Single Parameter
- `L`: Single Parameter

## Description
Cubic alpha `(α(T))` model. Also known as Twu-88 alpha.
```
αᵢ = Trᵢ^(N*(M-1))*exp(L*(1-Trᵢ^(N*M))
N = 2
Trᵢ = T/Tcᵢ
```
if `N` is specified, it will be used instead of the default value of 2.

## Model Construction Examples
```
# Using the default database
alpha = Twu88Alpha("water") #single input
alpha = Twu88Alpha(["water","ethanol"]) #multiple components

# Using user-provided parameters

# Passing files or folders
alpha = Twu88Alpha(["neon","hydrogen"]; userlocations = ["path/to/my/db","twu88.csv"])

# Passing parameters directly
alpha = Twu88Alpha(["neon","hydrogen"];
    userlocations = (;L = [0.40453, 156.21],
                    M = [0.95861, -0.0062072],
                    N = [0.8396, 5.047]) #if we don't pass N, then is assumed N = 2
                )
```

## References
1. Twu, C. H., Lee, L. L., & Starling, K. E. (1980). Improved analytical representation of argon thermodynamic behavior. Fluid Phase Equilibria, 4(1–2), 35–44. [doi:10.1016/0378-3812(80)80003-3](https://doi.org/10.1016/0378-3812(80)80003-3)
"""
function Twu88Alpha(components; userlocations = String[], verbose::Bool=false)
    _components = format_components(components)
    params = getparams(_components, ["alpha/Twu/Twu_like.csv"]; userlocations = userlocations, verbose = verbose,ignore_missing_singleparams = ["N"])
    M = params["M"]
    N = params["N"]
    L = params["L"]

    n = length(_components)
    for i in 1:n
        if N.ismissingvalues[i]
            N[i] = 2
        end
    end
    packagedparams = TwuAlphaParam(M,N,L)
    model = TwuAlpha(_components,packagedparams,String["10.1016/0378-3812(80)80003-3"])
    return model
end

function α_function(model::CubicModel,V,T,z,alpha_model::TwuAlphaModel)
    Tc = model.params.Tc.values
    _M  = alpha_model.params.M.values
    _N  = alpha_model.params.N.values
    _L  = alpha_model.params.L.values
    α = zeros(typeof(T*1.0),length(Tc))
    for i in @comps
        M = _M[i]
        N = _N[i]
        L = _L[i]
        Tr = T/Tc[i]
        α[i] = Tr^(N*(M-1))*exp(L*(1-Tr^(N*M)))
    end
    return α
end

function α_function(model::CubicModel,V,T,z::SingleComp,alpha_model::TwuAlphaModel)
    Tc = model.params.Tc.values[1]
    M  = alpha_model.params.M.values[1]
    N  = alpha_model.params.N.values[1]
    L  = alpha_model.params.L.values[1]
    Tr = T/Tc
    α = Tr^(N*(M-1))*exp(L*(1-Tr^(N*M)))
end

export TwuAlpha, Twu88Alpha, Twu91Alpha