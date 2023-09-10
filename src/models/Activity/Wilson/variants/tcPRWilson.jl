struct tcPRWilsonParam <: EoSParam
    g::PairParam{Float64}
    V::SingleParam{Float64}
end

struct tcPRWilson{c<:EoSModel} <: WilsonModel
    components::Array{String,1}
    params::tcPRWilsonParam
    puremodel::EoSVectorParam{c}
    references::Array{String,1}
end

export tcPRWilson

"""
    tcPRWilson <: WilsonModel
    tcPRWilson(components::Vector{String};
    puremodel = PR,
    userlocations = String[],
    pure_userlocations = String[],
    verbose = false)

## Input parameters
- `g`: Pair Parameter (`Float64`, asymetrical, defaults to `0`) - Interaction Parameter
- `v`: Single Parameter (`Float64`) - individual volumes.

## Input models
- `puremodel`: model to calculate pure pressure-dependent properties

## Description
tc-PR-Wilson activity model, meant to used in combination with the tc-PR cubic EoS:
```
Gᴱ = nRT∑xᵢlog(∑xⱼjΛᵢⱼ)
Λᵢⱼ = exp(-gᵢⱼ/T)*Vⱼ/Vᵢ
```
## References
1. Wilson, G. M. (1964). Vapor-liquid equilibrium. XI. A new expression for the excess free energy of mixing. Journal of the American Chemical Society, 86(2), 127–130. [doi:10.1021/ja01056a002](https://doi.org/10.1021/ja01056a002)
2. Piña-Martinez, A., Privat, R., Nikolaidis, I. K., Economou, I. G., & Jaubert, J.-N. (2021). What is the optimal activity coefficient model to be combined with the translated–consistent Peng–Robinson equation of state through advanced mixing rules? Cross-comparison and grading of the Wilson, UNIQUAC, and NRTL aE models against a benchmark database involving 200 binary systems. Industrial & Engineering Chemistry Research, 60(47), 17228–17247. [doi:10.1021/acs.iecr.1c03003](https://doi.org/10.1021/acs.iecr.1c03003)
"""
tcPRWilson

default_locations(::Type{tcPRWilson}) = ["Activity/tcPRWilson/tcPRWilson_unlike.csv"]


function tcPRWilson(components::Vector{String};
    puremodel = BasicIdeal,
    userlocations = String[],
    pure_userlocations = String[],
    verbose = false)

    params = getparams(components, default_locations(tcPRWilson); userlocations=userlocations, asymmetricparams=["g"], ignore_missing_singleparams=["g","V"], verbose=verbose)
    g = params["g"]
    V = params["V"]
    _puremodel = init_puremodel(puremodel,components,pure_userlocations,verbose)
    packagedparams = tcPRWilsonParam(g,V)
    references = String["10.1021/ja01056a002","10.1021/acs.iecr.1c03003"]
    model = tcPRWilson(components,packagedparams,_puremodel,references)
    return model
end

function wilson_volume(model::tcPRWilson,T)
    return model.params.V.values
end