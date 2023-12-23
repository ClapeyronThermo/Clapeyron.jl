abstract type CPCAMixingModel <: vdW1fRuleModel end

struct CPCAMixingParam <:EoSParam
    kb1::SingleParam{Float64}
    kb2::SingleParam{Float64}
    segment::SingleParam{Float64}
end

struct CPCAMixing <: CPCAMixingModel
    components::Array{String,1}
    params::CPCAMixingParam
    references::Array{String,1}
end

function CPCAMixing(components; activity = nothing, userlocations=String[],activity_userlocations=String[], verbose::Bool=false)
    params = getparams(components,String[], userlocations=userlocations, verbose=verbose)
    k1 = params["kb1"]
    k2 = params["kb2"]
    segment = params["segment"]
    param = CPCAMixingParam(kb1,kb2,segment)
    return CPCAMixing(format_components(components), param, ["10.1021/acs.iecr.2c00902"])
end

function mixing_rule(model::CPCAMixingModel, V, T, z, model_params)
    kb1 = model_params.kb1.values
    kb2 = model_params.kb2.values
    segment = model_params.segment.values
end

const CPCACubic = RK{I,CPAAlpha,T,CPCAMixing} where {I,T}

function cubic_ab(model::CPCACubic,V,T,z=SA[1.0],n=sum(z))
    a = model.params.a.values
    b = model.params.b.values
    T = T*float(one(T))
    α = @f(α_function,model.alpha)
    c = @f(translation,model.translation)
    if length(z) > 1
    ā,b̄,c̄ = @f(mixing_rule,model.mixing,α,a,b,c)
    else
        ā = a[1,1]*α[1]
        kb1 = model.mixing.params.kb1.values[1]
        kb2 = model.mixing.params.kb2.values[1]
        m = model.mixing.params.segment.values[1]
        Tc = model.params.Tc.values[1]
        β = kb1*exp(-kb2*T/Tc)
        b̄ = b[1,1]*β
        c̄ = c[1]
    end
    return ā ,b̄, c̄
end


abstract type CPCAModel <: CPAModel end

struct CPCAParam <: EoSParam
    Mw::SingleParam{Float64}
    Tc::SingleParam{Float64}
    segment::SingleParam{Float64}
    a::PairParam{Float64}
    b::PairParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64} 
end

struct CPCA{T <: IdealModel,c <: CubicModel} <: CPAModel
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
    CPCAModel <: CPAModel

    function CPCA(components;
        idealmodel=BasicIdeal,
        translation = NoTranslation
        userlocations=String[],
        ideal_userlocations=String[],
        translation_userlocations=String[],
        verbose=false,
        assoc_options = AssocOptions())

## Input parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `a`: Single Parameter (`Float64`) - Atraction parameter `[m^6*Pa/mol]`
- `b`: Single Parameter (`Float64`) - Covolume `[m^3/mol]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Paramater (no units)
- `l`: Pair Parameter (`Float64`) (optional) - Binary Interaction Paramater (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m^3]`

## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `a`: Pair Parameter (`Float64`) - Mixed Atraction Parameter `[m^6*Pa/mol]`
- `b`: Pair Parameter (`Float64`) - Mixed Covolume `[m^3/mol]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[J]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m^3]`

## Input models
- `idealmodel`: Ideal Model
- `translation`: Translation model

## Description

Cubic + Chain + Association (CPCA) EoS. Consists in the addition of a cubic part, a chain part and an association part:
```
a_res(model::CPA) = a_res(model::Cubic) + a_chain(model) + a_assoc(model)
```

## References
1. Kontogeorgis, G. M., Michelsen, M. L., Folas, G. K., Derawi, S., von Solms, N., & Stenby, E. H. (2006). Ten years with the CPA (cubic-plus-association) equation of state. Part 1. Pure compounds and self-associating systems. Industrial & Engineering Chemistry Research, 45(14), 4855–4868. [doi:10.1021/ie051305v](https://doi.org/10.1021/ie051305v)
"""
function CPCA(components;
    idealmodel=BasicIdeal,
    translation = NoTranslation
    userlocations=String[],
    ideal_userlocations=String[],
    translation_userlocations=String[],
    verbose=false,
    assoc_options = AssocOptions())
    locs = ["SAFT/CPA/CPCA/", "properties/molarmass.csv","properties/critical.csv"]
    params = getparams(components, locs; userlocations=userlocations, verbose=verbose)
    
    sites = get!(params,"sites") do
        SiteParam(components)
    end

    Pc = get!(params,"Pc") do
        SingleParam("Pc",components)
    end

    Mw  = params["Mw"]
    alpha
    
    
    k = get(params,"k",nothing)
    Tc = params["Tc"]
    a  = epsilon_LorentzBerthelot(params["a"], k)
    b  = sigma_LorentzBerthelot(params["b"], l)

    epsilon_assoc = get!(params,"epsilon_assoc") do
        AssocParam("epsilon_assoc",components)
    end

    bondvol = get!(params,"bondvol") do
        AssocParam("bondvol",components)
    end

    bondvol,epsilon_assoc = assoc_mix(bondvol,epsilon_assoc,cbrt.(b),assoc_options)
    packagedparams = CPCAParam(a, b, Tc, epsilon_assoc, bondvol, Mw)
    
    #init cubic model
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    init_alpha = CPAAlpha(components,userlocations = (;c1 = 10.2425 .* Mw .^ 0.2829)) 
    init_mixing = CPCAMixing(components,userlocations = (;c1 = 0.0, c2 = 0.0))
    init_translation = init_model(translation,components,translation_userlocations,verbose)
    cubicparams = ABCubicParam(a, b, params["Tc"],Pc,Mw) #PR, RK, vdW
    init_cubicmodel = cubicmodel(components,init_alpha,init_mixing,init_translation,cubicparams,init_idealmodel,String[])

    references = ["10.1021/ie051305v"]
    model = CPCA(components, radial_dist, init_cubicmodel, packagedparams, sites, init_idealmodel, assoc_options, references)
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

function a_res(model::CPAModel, V, T, z, _data = @f(data))
    n,ā,b̄,c̄ = _data
    return a_res(model.cubicmodel,V,T,z,_data) + a_assoc(model,V+c̄*n,T,z,_data)
end
