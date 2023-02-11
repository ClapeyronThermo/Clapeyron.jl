abstract type sCPAModel <: CPAModel end

struct sCPA{T <: IdealModel,c <: CubicModel} <: sCPAModel
    components::Array{String,1}
    cubicmodel::c
    params::CPAParam
    sites::SiteParam
    idealmodel::T
    assoc_options::AssocOptions
    references::Array{String,1}
end

@registermodel sCPA
export sCPA

"""
    sCPAModel <: CPAModel

    function sCPA(components; 
        idealmodel=BasicIdeal, 
        cubicmodel=RK, 
        alpha=sCPAAlpha, 
        mixing=vdW1fRule,
        activity=nothing,
        translation=NoTranslation, 
        userlocations=String[], 
        ideal_userlocations=String[], 
        alpha_userlocations=String[],
        activity_userlocations=String[],
        mixing_userlocations=String[],
        translation_userlocations=String[],
        verbose=false,
        assoc_options = AssocOptions())

## Input parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `a`: Single Parameter (`Float64`) - Atraction Parameter
- `b`: Single Parameter (`Float64`) - Covolume
- `c1`: Single Parameter (`Float64`) - α-function constant Parameter
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Paramater (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m^3]`

## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `a`: Pair Parameter (`Float64`) - Mixed Atraction Parameter
- `b`: Pair Parameter (`Float64`) - Mixed Covolume
- `c1`: Single Parameter (`Float64`) - α-function constant Parameter
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m^3]`

## Input models
- `idealmodel`: Ideal Model
- `cubicmodel`: Cubic Model

## Description

simplified CPA

## References
1. Kontogeorgis, G. M., Michelsen, M. L., Folas, G. K., Derawi, S., von Solms, N., & Stenby, E. H. (2006). Ten years with the CPA (cubic-plus-association) equation of state. Part 1. Pure compounds and self-associating systems. Industrial & Engineering Chemistry Research, 45(14), 4855–4868. [doi:10.1021/ie051305v](https://doi.org/10.1021/ie051305v)
"""
sCPA

function sCPA(components; 
            idealmodel=BasicIdeal, 
            cubicmodel=RK, 
            alpha=sCPAAlpha, 
            mixing=vdW1fRule,
            activity=nothing,
            translation=NoTranslation, 
            userlocations=String[], 
            ideal_userlocations=String[], 
            alpha_userlocations=String[],
            activity_userlocations=String[],
            mixing_userlocations=String[],
            translation_userlocations=String[],
            verbose=false,
            assoc_options = AssocOptions(), kwargs...)

    params,sites = getparams(components, ["SAFT/CPA/sCPA/", "properties/molarmass.csv","properties/critical.csv"]; userlocations=userlocations, verbose=verbose)
    Mw  = params["Mw"]
    k = get(params,"k",nothing)
    Tc = params["Tc"]
    c1 = params["c1"]
    params["a"].values .*= 1E-1
    params["b"].values .*= 1E-3
    a  = epsilon_LorentzBerthelot(params["a"], k)
    b  = sigma_LorentzBerthelot(params["b"])
    epsilon_assoc = params["epsilon_assoc"]
    bondvol = params["bondvol"]
    bondvol,epsilon_assoc = assoc_mix(bondvol,epsilon_assoc,cbrt.(b),assoc_options)
    packagedparams = CPAParam(a, b, c1, Tc, epsilon_assoc, bondvol,Mw)

    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    init_alpha = init_model(alpha,components,alpha_userlocations,verbose)
    init_mixing = init_model(mixing,components,activity,mixing_userlocations,activity_userlocations,verbose)
    init_translation = init_model(translation,components,translation_userlocations,verbose)

    if occursin("RK",string(cubicmodel))
        cubicparams = RKParam(a, b, params["Tc"],params["Pc"],Mw)
    elseif occursin("PR",string(cubicmodel))
        cubicparams = PRParam(a, b, params["Tc"],params["Pc"],Mw)
    end

    init_cubicmodel = cubicmodel(components,init_alpha,init_mixing,init_translation,cubicparams,init_idealmodel,String[])

    references = ["10.1021/ie051305v"]

    model = sCPA(components, init_cubicmodel, packagedparams, sites, init_idealmodel, assoc_options, references)
    return model
end

function Δ(model::sCPAModel, V, T, z, i, j, a, b)
    ϵ_associjab = model.params.epsilon_assoc.values[i,j][a,b] * 1e2/R̄
    βijab = model.params.bondvol.values[i,j][a,b] * 1e-3
    Σz = ∑(z)
    b = model.params.b.values
    b̄ = dot(z,b,z)/(Σz*Σz)
    η = b̄*Σz/(4*V)
    g = (1-1.9η)^-1
    bij = (b[i,i]+b[j,j])/2
    return g*(exp(ϵ_associjab/T)-1)*βijab*bij/N_A
end
