struct CPAParam <: EoSParam
    Mw::SingleParam{Float64}
    Tc::SingleParam{Float64}
    a::PairParam{Float64}
    b::PairParam{Float64}
    c1::SingleParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
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

function transform_params(::Type{CPA},params,_components)
    a = PairParam(params["a"])
    _components = a.components
    b = PairParam(params["b"])
    params["a"] = a
    params["b"] = b
    
    sites = get!(params,"sites") do
        SiteParam(_components)
    end

    Mw  = params["Mw"]
    Tc = params["Tc"]
    c1 = get(params,"c1",nothing)
    Pc = get!(params,"Pc") do
        SingleParam("Pc",_components)
    end
    params["Pc"] = Pc
    epsilon_assoc = get!(params,"epsilon_assoc") do
        AssocParam("epsilon_assoc",_components)
    end
    bondvol = get!(params,"bondvol") do
        AssocParam("bondvol",_components)
    end

    #assoc mixing is done later.
    #bondvol,epsilon_assoc = assoc_mix(bondvol,epsilon_assoc,cbrt.(b),assoc_options,sites)
    params
end

"""
    CPAModel <: EoSModel

    function CPA(components;
        radial_dist::Symbol = :CS,
        idealmodel = BasicIdeal,
        cubicmodel = RK,
        alpha = sCPAAlpha,
        mixing = vdW1fRule,
        activity = nothing,
        translation = NoTranslation,
        userlocations = String[],
        ideal_userlocations = String[],
        alpha_userlocations = String[],
        activity_userlocations = String[],
        mixing_userlocations = String[],
        translation_userlocations = String[],
        reference_state = nothing,
        verbose = false,
        assoc_options = AssocOptions())

## Input parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `a`: Single Parameter (`Float64`) - Atraction parameter `[m^6*Pa/mol]`
- `b`: Single Parameter (`Float64`) - Covolume `[m^3/mol]`
- `c1`: Single Parameter (`Float64`) - α-function constant Parameter (no units)
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Paramater (no units)
- `l`: Pair Parameter (`Float64`) (optional) - Binary Interaction Paramater (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m^3]`

## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
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

default_references(::Type{CPA}) = ["10.1021/ie051305v"]
default_locations(::Type{CPA}) = ["SAFT/CPA", "properties/molarmass.csv","properties/critical.csv"]

function default_locations(model::CPA)
    locs = if radial_dist == :CS
        ["SAFT/CPA", "properties/molarmass.csv","properties/critical.csv"]
    elseif radial_dist == :KG
        ["SAFT/CPA/sCPA/", "properties/molarmass.csv","properties/critical.csv"]
    else
        throw(error("CPA: incorrect specification of radial_dist, try using `:CS` (original CPA) or `:KG` (simplified CPA)"))
    end
end

export CPA
function CPA(components;
    idealmodel = BasicIdeal,
    radial_dist::Symbol = :CS,
    cubicmodel = RK,
    alpha = CPAAlpha,
    mixing = vdW1fRule,
    activity = nothing,
    translation = NoTranslation,
    userlocations = String[],
    ideal_userlocations = String[],
    alpha_userlocations = String[],
    activity_userlocations = String[],
    mixing_userlocations = String[],
    translation_userlocations = String[],
    reference_state = nothing,
    verbose = false,
    assoc_options = AssocOptions())

    locs = if radial_dist == :CS
        ["SAFT/CPA", "properties/molarmass.csv","properties/critical.csv"]
    elseif radial_dist == :KG
        ["SAFT/CPA/sCPA/", "properties/molarmass.csv","properties/critical.csv"]
    else
        throw(error("CPA: incorrect specification of radial_dist, try using `:CS` (original CPA) or `:KG` (simplified CPA)"))
    end

    _components = format_components(components)
    params = getparams(_components, locs; userlocations = userlocations, verbose = verbose, ignore_missing_singleparams = ["Pc","Vc","acentricfactor"])    
    references = ["10.1021/ie051305v"]
    
    #adds missing parameters, in the correct shape
    transform_params(CPA,params,components)
    
    #creates a CPAParam
    packagedparams = CPAParam(params["Mw"],params["Tc"],params["a"],params["b"],params["c1"],params["epsilon_assoc"],params["bondvol"])
    
    #init empty cubic model
    init_cubicmodel = CubicModel(cubicmodel,params,components;
                                idealmodel,alpha,mixing,activity,translation,
                                userlocations,ideal_userlocations,alpha_userlocations,activity_userlocations,mixing_userlocations,translation_userlocations,
                                reference_state, verbose)

    init_idealmodel = init_cubicmodel.idealmodel

    model = CPA(_components, radial_dist, init_cubicmodel, packagedparams, params["sites"], init_idealmodel, assoc_options, references)
    
    #perform mixing rules for a and b, perform association combining rules.
    k = get(params,"k",nothing)
    l = get(params,"l",nothing)
    recombine_cpa!(model,k,l)

    #calculates reference states, if any.
    set_reference_state!(model,reference_state;verbose)
    return model
end

#needs to add overload for other mixing rules
function ab_premixing(model::CPAModel,mixing::MixingRule, k, l)
    a = model.params.a
    b = model.params.b
    epsilon_LorentzBerthelot!(a,k)
    sigma_LorentzBerthelot!(b,l)
    return a,b
end


function recombine_cpa!(model::CPAModel,k = nothing, l = nothing)
    a = model.params.a
    b = model.params.b
    cubicmodel = model.cubicmodel
    Pc = cubicmodel.params.Pc
    components = model.components
    for i in 1:length(components)
        if Pc.ismissingvalues[i]
            Ωa,Ωb = ab_consts(cubicmodel)
            Ωa,Ωb = ab_consts(cubicmodel)
            b̄r = b[i,i]/Ωb
            ār = a[i,i]/Ωa
            Pc[i] = ār/(b̄r*b̄r)
        end
    end

    #recombine_mixing also calculates additional things for some ABC cubics
    recombine_mixing!(model,cubicmodel.mixing,k,l)
    cubicmodel.params.a.values .= a.values
    cubicmodel.params.b.values .= b.values
    recombine_translation!(cubicmodel,cubicmodel.translation)
    
    if model.cubicmodel.alpha isa CPAAlphaModel
        cubicmodel.alpha.params.c1.values .= model.params.c1.values
    else
        recombine_alpha!(cubicmodel,cubicmodel.alpha)
    end

    #association recombining
    recombine_assoc!(model,cbrt.(b))

    return model
end

function recombine_impl!(model::CPAModel)
    #this function recombines all parameters of CPA.
    #for polar CPA or other variants, more code needs to be added here.
    recombine_cpa!(model)

    return model
end

lb_volume(model::CPAModel,T,z) = lb_volume(model.cubicmodel,T,z)
T_scale(model::CPAModel,z) = T_scale(model.cubicmodel,z)

ab_consts(model::CPAModel) = ab_consts(model.cubicmodel)
ab_consts(model::CPAModel,z) = ab_consts(model.cubicmodel,z)
ab_consts(::Type{T}) where T <: CPAModel = ab_consts(fieldtype(T,:cubicmodel))

function p_scale(model::CPAModel,z)
    #does not depend on Pc, so it can be made optional on CPA input
    b = model.cubicmodel.params.b.values
    a = model.cubicmodel.params.a.values
    Ωa,Ωb = ab_consts(model,z)
    b̄r = dot(z,b,z)/(sum(z)*Ωb)
    ār = dot(z,a,z)/Ωa
    return ār/(b̄r*b̄r)
end

function cpa_is_pure_cubic(model::CPAModel)
    assoc_pair_length(model) == 0
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
    z = SA[1.0]
    T = T_scale(model,z)
    lb_v = lb_volume(model,T,z)
    return (1.0, log10(lb_v/0.3))
end

function crit_pure(model::CPAModel)
    if cpa_is_pure_cubic(model) && !model.cubicmodel.params.Pc.ismissingvalues[1]
        return crit_pure(model.cubicmodel)
    else
        return crit_pure(model,x0_crit_pure(model))
    end
end

function x0_sat_pure(model::CPAModel,T,crit = nothing)
    #use the cubic initial guess if we don't have association.
    cpa_is_pure_cubic(model) && x0_sat_pure(model.cubicmodel,T)
    if crit === nothing
        _,vl,vv = x0_sat_pure_virial(model,T)
    else
        _,vl,vv = x0_sat_pure_crit(model,T,crit)
    end
    return vl,vv
end

function x0_psat(model::CPAModel,T,crit = nothing)
    cpa_is_pure_cubic(model) && x0_sat_pure(model.cubicmodel,T)
    if crit === nothing
        p,_,_ = x0_sat_pure_virial(model,T)
    else
        p,_,_ = x0_sat_pure_crit(model,T,crit)
    end
    return p
end

#=
if we don't have association, reduce to the inner cubic model.
=#
function volume_impl(model::CPAModel,p,T,z,phase,threaded,vol0)
    if cpa_is_pure_cubic(model)
        return volume_impl(model.cubicmodel,p,T,z,phase,threaded,vol0)
    else
        return default_volume_impl(model,p,T,z,phase,threaded,vol0)
    end
end

#approximating the gas phase as the pure cubic.
function x0_volume_gas(model::CPAModel, p, T, z)
    return volume(model.cubicmodel,p,T,z,phase = :v)
end

function x0_volume_liquid(model::CPAModel,p, T, z)
    cpa_is_pure_cubic(model) && return volume(model.cubicmodel,p,T,z,phase = :l)
    return 1.1*lb_volume(model,T,z)
end

data(model::CPAModel, V, T, z) = data(model.cubicmodel,V,T,z)

function a_res(model::CPAModel, V, T, z, _data = @f(data))
    n,ā,b̄,c̄ = _data
    return a_res(model.cubicmodel,V,T,z,_data) + a_assoc(model,V0 + c̄*n,T,z,_data)
end

function Δ(model::CPAModel, V, T, z, i, j, a, b, _data = @f(data))
    n,ā,b̄,c̄ = _data
    ϵ_associjab = model.params.epsilon_assoc.values[i,j][a,b]
    βijab = model.params.bondvol.values[i,j][a,b]
    b = model.params.b.values
    η = n*b̄/(4*V)
    rdf = model.radial_dist
    if rdf == :CS #CPA original
        g = (1-0.5*η)/(1-η)^3
    elseif rdf == :KG #sCPA
        g = 1/(1-1.9η)
    else
        g = zero(η)/zero(η)
    end

    return g*expm1(ϵ_associjab/T)*βijab*b[i,j]/N_A
end

#optimized Δ function for CPA, we only calculate g once.
function  Δ(model::CPAModel, V, T, z,_data=@f(data))
    n,ā,b̄,c̄ = _data
    β = model.params.bondvol.values
    b_cubic = model.params.b.values
    η = n*b̄/(4*V)
    rdf = model.radial_dist
    if rdf == :CS #CPA original
        g = (1-0.5*η)/(1-η)^3
    elseif rdf == :KG #sCPA
        g = 1/(1-1.9η)
    else
        g = zero(η)/zero(η)
    end
    Δout = assoc_similar(β,typeof(V+T+first(z)+one(eltype(model))))
    ϵ_assoc = model.params.epsilon_assoc
    Δout.values .= false  #fill with zeros, maybe it is not necessary?
    for (idx,(i,j),(a,b)) in indices(Δout)
        βijab = β[idx]
        if βijab != 0
            Δout[idx] = g*expm1(ϵ_assoc[i,j][a,b]/T)*βijab*b_cubic[i,j]/N_A
        end
    end
    return Δout
end
