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
    params = getparams(formatted_components, default_locations(Wilson); userlocations = userlocations, asymmetricparams=["g"], ignore_missing_singleparams=["g"], verbose = verbose)
    A₁₂        = params["A₁₂"]
    A₂₁        = params["A₂₁"]
    Mw        = params["Mw"]
    
    _puremodel = init_puremodel(puremodel,components,pure_userlocations,verbose)
    packagedparams = MargulesParam(A₁₂,A₂₁,Mw)
    references = String["10.1021/ja01056a002"]
    model = Margules(formatted_components,packagedparams,_puremodel,references)
    set_reference_state!(model,reference_state,verbose = verbose)
    binary_component_check(method,model)
    return model
end


#=
function activity_coefficient(model::MargulesModel,p,T,z)
    return activity_coefficient_margules(model,p,T,z)
end

function activity_coefficient_wilson(model::MargulesModel,p,T,z)
    Λ = (Vi' ./ Vi) .*exp.(-model.params.g.values/R̄/T)
    x = z ./ sum(z)
    lnγ = 1 .- log.(sum(x[i]*Λ[:,i] for i ∈ @comps)) .-sum(x[j] .*Λ[j,:] ./(sum(x[i]*Λ[j,i] for i ∈ @comps)) for j ∈ @comps)
    return exp.(lnγ)
end =#


function excess_g_wilson(model::WilsonModel,p,T,z,V = wilson_volume(model,T))
    g = model.params.g.values
    _0 = zero(T+first(z))
    n = sum(z)
    invn = 1/n
    invRT = 1/(R̄*T)
    res = _0
    for i ∈ @comps
        ∑xΛ = _0
        xi = z[i]*invn
        for j ∈ @comps
            Λij = exp(-g[i,j]*invRT)*V[j]/V[i]

            ∑xΛ += Λij*z[j]*invn
        end
        res += xi*log(∑xΛ)
    end
    return -n*res*R̄*T
end