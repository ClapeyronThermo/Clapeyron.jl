struct CPAParam <: EoSParam
    a::PairParam{Float64}
    b::PairParam{Float64}
    c1::SingleParam{Float64}
    Tc::SingleParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
    Mw::SingleParam{Float64}
end


struct CPA{T <: IdealModel,c <: CubicModel} <: CPAModel
    components::Array{String,1}
    radial_dist::Symbol
    cubicmodel::c
    params::CPAParam
    sites::SiteParam
    idealmodel::T
    assoc_options::AssocOptions
    references::Array{String,1}
end

"""
    CPAModel <: EoSModel

    function CPA(components;
        radial_dist::Symbol = :CS,
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
- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Pressure `[Pa]`
- `a`: Single Parameter (`Float64`) - Atraction parameter `[m^6*Pa/mol]`
- `b`: Single Parameter (`Float64`) - Covolume `[m^3/mol]`
- `c1`: Single Parameter (`Float64`) - α-function constant Parameter (no units)
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Paramater (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m^3]`

## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `a`: Pair Parameter (`Float64`) - Mixed Atraction Parameter `[m^6*Pa/mol]`
- `b`: Pair Parameter (`Float64`) - Mixed Covolume `[m^3/mol]`
- `c1`: Single Parameter (`Float64`) - α-function constant Parameter (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[J]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m^3]`

## Input models
- `idealmodel`: Ideal Model
- `cubicmodel`: Cubic Model

## Description

Cubic Plus Association (CPA) EoS. Consists in the addition of a cubic part and an association part:
```
a_res(model::CPA) = a_res(model::Cubic) + a_assoc(model)
```

The `radial_dist` argument can be used to choose between a Carnahan-Starling form (`CS`, default) or the Kontogeorgis (`KG`) term, more widely known as s-CPA.

## References
1. Kontogeorgis, G. M., Michelsen, M. L., Folas, G. K., Derawi, S., von Solms, N., & Stenby, E. H. (2006). Ten years with the CPA (cubic-plus-association) equation of state. Part 1. Pure compounds and self-associating systems. Industrial & Engineering Chemistry Research, 45(14), 4855–4868. [doi:10.1021/ie051305v](https://doi.org/10.1021/ie051305v)
"""
CPA

export CPA
function CPA(components;
    idealmodel=BasicIdeal,
    radial_dist::Symbol = :CS,
    cubicmodel=RK,
    alpha=CPAAlpha,
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

    locs = if radial_dist == :CS
        ["SAFT/CPA", "properties/molarmass.csv","properties/critical.csv"]
    elseif radial_dist == :KG
        ["SAFT/CPA/sCPA/", "properties/molarmass.csv","properties/critical.csv"]
    else
        throw(error("CPA: incorrect specification of radial_dist, try using `:CS` (original CPA) or `:KG` (simplified CPA)"))
    end

    params = getparams(components, locs; userlocations=userlocations, verbose=verbose)
    
    sites = get!(params,"sites") do
        SiteParam(components)
    end

    Pc = get!(params,"Pc") do
        SingleParam("Pc",components)
    end

    Mw  = params["Mw"]
    k = get(params,"k",nothing)
    Tc = params["Tc"]
    c1 = params["c1"]
    a  = epsilon_LorentzBerthelot(params["a"], k)
    b  = sigma_LorentzBerthelot(params["b"])

    epsilon_assoc = get!(params,"epsilon_assoc") do
        AssocParam("epsilon_assoc",components)
    end

    bondvol = get!(params,"bondvol") do
        AssocParam("bondvol",components)
    end

    bondvol,epsilon_assoc = assoc_mix(bondvol,epsilon_assoc,cbrt.(b),assoc_options)
    packagedparams = CPAParam(a, b, c1, Tc, epsilon_assoc, bondvol, Mw)
    
    #init cubic model
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    init_alpha = init_model(alpha,components,alpha_userlocations,verbose)
    init_mixing = init_model(mixing,components,activity,mixing_userlocations,activity_userlocations,verbose)
    init_translation = init_model(translation,components,translation_userlocations,verbose)
    cubicparams = ABCubicParam(a, b, params["Tc"],Pc,Mw) #PR, RK, vdW
    init_cubicmodel = cubicmodel(components,init_alpha,init_mixing,init_translation,cubicparams,init_idealmodel,String[])

    references = ["10.1021/ie051305v"]

    model = CPA(components, radial_dist, init_cubicmodel, packagedparams, sites, init_idealmodel, assoc_options, references)
    return model
end

function recombine_impl!(model::CPAModel)
    assoc_options = model.assoc_options
    a = model.params.a
    b = model.params.b

    a  = epsilon_LorentzBerthelot!(a)
    b  = sigma_LorentzBerthelot!(b)

    epsilon_assoc = model.params.epsilon_assoc
    bondvol = model.params.bondvol
    bondvol,epsilon_assoc = assoc_mix(bondvol,epsilon_assoc,cbrt.(b),assoc_options) #combining rules for association

    model.params.epsilon_assoc.values.values[:] = epsilon_assoc.values.values
    model.params.bondvol.values.values[:] = bondvol.values.values
    return model
end

lb_volume(model::CPAModel,z = SA[1.0]) = lb_volume(model.cubicmodel,z)
T_scale(model::CPAModel,z=SA[1.0]) = T_scale(model.cubicmodel,z)
function p_scale(model::CPAModel,z=SA[1.0])
    #does not depend on Pc, so it can be made optional on CPA input
    b = model.cubicmodel.params.b.values
    a = model.cubicmodel.params.a.values
    Ωa,Ωb = ab_consts(model.cubicmodel)
    b̄r = dot(z,b,z)/(sum(z)*Ωb)
    ār = dot(z,a,z)/Ωa
    return ār/(b̄r*b̄r)
end

function show_info(io,model::CPAModel) 
    rdf = model.radial_dist
    println(io)
    if rdf == :CS #CPA original
        print(io,"RDF: Carnahan-Starling (original CPA)")
    elseif rdf == :KG || rdf == :OT #sCPA
        print(io,"RDF: Kontogeorgis (s-CPA)")
    end
end

function x0_crit_pure(model::CPAModel)
    lb_v = lb_volume(model)
    return (1.0, log10(lb_v/0.3))
end

data(model::CPAModel, V, T, z) = data(model.cubicmodel,V,T,z)

function a_res(model::CPAModel, V, T, z, _data = @f(_data))
    n,ā,b̄,c̄ = _data
    return a_res(model.cubicmodel,V,T,z,_data) + a_assoc(model,V+c̄*n,T,z,_data)
end

ab_consts(model::CPAModel) = ab_consts(model.cubicmodel)

function Δ(model::CPAModel, V, T, z, i, j, a, b, _data = data(model.cubicmodel,V,T,z))
    n,ā,b̄,c̄ = _data
    ϵ_associjab = model.params.epsilon_assoc.values[i,j][a,b]
    βijab = model.params.bondvol.values[i,j][a,b]
    b = model.params.b.values
    η = n*b̄/(4*V)

    rdf = model.radial_dist
    g = if rdf == :CS #CPA original
        (1-0.5*η)/(1-η)^3
    elseif rdf == :KG #sCPA
        1/(1-1.9η)
    else
        zero(η)/zero(η)
    end

    return g*expm1(ϵ_associjab/T)*βijab*b[i,j]/N_A
end
