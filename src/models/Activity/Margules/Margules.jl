struct MargulesParam <: EoSParam
    A₁₂::SingleParam{Float64}
    A₂₁::SingleParam{Float64}
    Mw::SingleParam{Float64}
end

abstract type MargulesModel <: ActivityModel end

struct Margules{c<:EoSModel} <: MargulesModel
    components::Array{String,1}
    params::MargulesParam
    puremodel::EoSVectorParam{c}
    references::Array{String,1}
end

export Margules

"""
    Margules <: ActivityModel
    Margules(components;
    puremodel = PR,
    userlocations = String[],
    pure_userlocations = String[],
    verbose = false,
    reference_state = nothing)

## Input parameters
- `A₁₂`: Single Parameter (`Float64`, defaults to `0`) - Binary Interaction Parameter
- `A₂₁`: Single Parameter (`Float64`, defaults to `0`) - Binary Interaction Parameter
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`

## model parameters
- `A₁₂`: Single Parameter (`Float64`, defaults to `0`) - Binary Interaction Parameter
- `A₂₁`: Single Parameter (`Float64`, defaults to `0`) - Binary Interaction Parameter

## Input models
- `puremodel`: model to calculate pure pressure-dependent properties

## Description
Margules activity coefficient model, for binary mixture:
```
Gᴱ = nRT·(x₁·x₂·(A₂₁·x₁+A₁₂·x₂))
```

## Model Construction Examples
```
# Using the default database
model = Margules(["water","ethanol"]) #Default pure model: PR
model = Margules(["water","ethanol"],puremodel = BasicIdeal) #Using Ideal Gas for pure model properties
model = Margules(["water","ethanol"],puremodel = PCSAFT) #Using Real Gas model for pure model properties

# Passing a prebuilt model

my_puremodel = AbbottVirial(["water","ethanol"]; userlocations = ["path/to/my/db","critical.csv"])
mixing = Margules(["water","ethanol"],puremodel = my_puremodel)

# Using user-provided parameters

# Passing files or folders
model = Margules(["water","ethanol"];userlocations = ["path/to/my/db","margules.csv"])

# Passing parameters directly
model = Margules(["water","ethanol"],
        userlocations = (A₁₂ = [4512],
                        A₂₁ = [3988.52],
                        Mw = [18.015, 46.069])
                        )
```

## References
1. Max Margules, « Über die Zusammensetzung der gesättigten Dämpfe von Misschungen », Sitzungsberichte der Kaiserliche Akadamie der Wissenschaften Wien Mathematisch-Naturwissenschaftliche Klasse II, vol. 104, 1895, p. 1243–1278
"""
Margules

default_locations(::Type{Margules}) = ["properties/critical.csv", "properties/molarmass.csv","Activity/Margules/margules_unlike.csv"]

function Margules(components;
    puremodel = PR,
    userlocations = String[],
    pure_userlocations = String[],
    verbose = false,
    reference_state = nothing)
    
    formatted_components = format_components(components)
    params = getparams(formatted_components, default_locations(Margules); userlocations = userlocations, asymmetricparams=["A₁₂","A₂₁"], ignore_missing_singleparams=["A₁₂","A₂₁","Mw"], verbose = verbose)
    A₁₂        = params["A₁₂"]
    A₂₁        = params["A₂₁"]
    Mw         = params["Mw"]
    
    _puremodel = init_puremodel(puremodel,components,pure_userlocations,verbose)
    packagedparams = MargulesParam(A₁₂,A₂₁,Mw)
    references = String[""]
    model = Margules(formatted_components,packagedparams,_puremodel,references)
    set_reference_state!(model,reference_state,verbose = verbose)
    binary_component_check(Margules,model)
    return model
end

function excess_g_margules(model::MargulesModel, p, T, z)
    n = sum(z)
    x = z ./ n

    A₁₂ = model.params.A₁₂.values[1]
    A₂₁ = model.params.A₂₁.values[1]

    ge = x[1] * x[2] * (A₂₁*x[1] + A₁₂*x[2])

    return n * R̄ * T * ge
end