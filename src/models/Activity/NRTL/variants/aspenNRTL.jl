struct aspenNRTLParam <: EoSParam
    a0::PairParam{Float64}
    a1::PairParam{Float64}
    t0::PairParam{Float64}
    t1::PairParam{Float64}
    t2::PairParam{Float64}
    t3::PairParam{Float64}
end

abstract type aspenNRTLModel <: NRTLModel end

struct aspenNRTL{c<:EoSModel} <: aspenNRTLModel
    components::Array{String,1}
    params::aspenNRTLParam
    puremodel::EoSVectorParam{c}
    references::Array{String,1}
end

export aspenNRTL
"""
    aspenNRTL <: ActivityModel

    function aspenNRTL(components;
    puremodel=PR,
    userlocations = String[],
    pure_userlocations = String[],
    verbose = false)

## Input parameters
- `a0`: Pair Parameter (`Float64`, asymetrical, defaults to `0`) - Interaction Parameter
- `a1`: Pair Parameter (`Float64`, asymetrical, defaults to `0`) - Interaction Parameter
- `t0`: Pair Parameter (`Float64`, asymetrical, defaults to `0`) - Interaction Parameter
- `t1`: Pair Parameter (`Float64`, asymetrical, defaults to `0`) - Interaction Parameter
- `t2`: Pair Parameter (`Float64`, asymetrical, defaults to `0`) - Interaction Parameter
- `t3`: Pair Parameter (`Float64`, asymetrical, defaults to `0`) - Interaction Parameter

## Input models
- `puremodel`: model to calculate pure pressure-dependent properties

## Description
NRTL (Non Random Two Fluid) activity model:
```
Gᴱ = nRT∑[xᵢ(∑τⱼᵢGⱼᵢxⱼ)/(∑Gⱼᵢxⱼ)]
Gᵢⱼ exp(-αᵢⱼτᵢⱼ)
αᵢⱼ = αᵢⱼ₀ + αᵢⱼ₁T
τᵢⱼ = tᵢⱼ₀ + tᵢⱼ₁/T + tᵢⱼ₂*ln(T) + tᵢⱼ₃*T
```

## Model Construction Examples
```
# Using the default database
model = aspenNRTL(["water","ethanol"]) #Default pure model: PR
model = aspenNRTL(["water","ethanol"],puremodel = BasicIdeal) #Using Ideal Gas for pure model properties
model = aspenNRTL(["water","ethanol"],puremodel = PCSAFT) #Using Real Gas model for pure model properties
# Passing a prebuilt model

my_puremodel = AbbottVirial(["water","ethanol"]; userlocations = ["path/to/my/db","critical.csv"])
mixing = aspenNRTL(["water","ethanol"],puremodel = my_puremodel)

# Using user-provided parameters

# Passing files or folders
model = aspenNRTL(["water","ethanol"];userlocations = ["path/to/my/db","nrtl.csv"])

# Passing parameters directly
model = aspenNRTL(["water","acetone"],
userlocations = (a0 = [0.0 0.228632; 0.228632 0.0],
                a1 = [0.0 0.000136; 0.000136 0.0],
                t0 = [0.0 13.374756; 16.081848 0.0],
                t1 = [0.0 -415.527344; -88.8125 0.0],
                t2 = [0.0 -1.913689; -2.930901 0.0],
                t3 = [0.0 0.00153; 0.005305 0.0])
            )
```

## References
1. Renon, H., & Prausnitz, J. M. (1968). Local compositions in thermodynamic excess functions for liquid mixtures. AIChE journal. American Institute of Chemical Engineers, 14(1), 135–144. [doi:10.1002/aic.690140124](https://doi.org/10.1002/aic.690140124)
"""
aspenNRTL

default_locations(::Type{aspenNRTL}) = ["Activity/NRTL/aspenNRTL/aspenNRTL_unlike.csv"]

function aspenNRTL(components; puremodel=PR,
    userlocations = String[], 
    pure_userlocations = String[],
    verbose = false)

    formatted_components = format_components(components)
    params = getparams(formatted_components, default_locations(aspenNRTL); userlocations = userlocations, asymmetricparams=["t0","t1","t2","t3"], ignore_missing_singleparams=asymmetricparams=["t0","t1","t2","t3"], verbose = verbose)
    a0  = params["a0"]
    a1  = params["a1"]
    t0  = params["t0"]
    t1  = params["t1"]
    t2  = params["t2"]
    t3  = params["t3"]

    _puremodel = init_puremodel(puremodel,components,pure_userlocations,verbose)
    packagedparams = aspenNRTLParam(a0,a1,t0,t1,t2,t3)
    references = String["10.1002/aic.690140124"]
    model = aspenNRTL(formatted_components,packagedparams,_puremodel,references)
    return model
end

function aspenNRTL(model::NRTL)
    params = model.params
    a0 = deepcopy(params.c)
    a1 = deepcopy(params.a)
    a1 .= 0
    t0 = deepcopy(params.a)
    t1 = deepcopy(params.b)
    t2 = deepcopy(params.a)
    t3 = deepcopy(params.a)
    t2 .= 0
    t3 .= 0
    packagedparams = aspenNRTLParam(a0,a1,t0,t1,t2,t3)
    return aspenNRTL(model.components,packagedparams,model.puremodel,model.references)
end
function excess_g_res(model::aspenNRTLModel,p,T,z)
    a₀ = model.params.a0.values
    a₁  = model.params.a1.values
    t₀  = model.params.t0.values
    t₁  = model.params.t1.values
    t₂  = model.params.t2.values
    t₃  = model.params.t3.values

    _0 = zero(T+first(z))
    n = sum(z)
    invn = 1/n
    invT = 1/(T)
    lnT = log(T)
    res = _0 
    for i ∈ @comps
        ∑τGx = _0
        ∑Gx = _0
        xi = z[i]*invn
        for j ∈ @comps
            xj = z[j]*invn
            α = a₀[j,i] + a₁[j,i]*T
            τji = t₀[j,i] + t₁[j,i]*invT + t₂[j,i]*lnT + t₃[j,i]*T
            Gji = exp(-α*τji)
            Gx = xj*Gji
            ∑Gx += Gx
            ∑τGx += Gx*τji
        end
        res += xi*∑τGx/∑Gx
    end
    return n*res*R̄*T
end