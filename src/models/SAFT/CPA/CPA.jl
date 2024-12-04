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
    params = getparams(_components, locs; userlocations = userlocations, verbose = verbose)
    
    sites = get!(params,"sites") do
        SiteParam(_components)
    end

    Pc = get!(params,"Pc") do
        SingleParam("Pc",_components)
    end

    Mw  = params["Mw"]
    k = get(params,"k",nothing)
    l = get(params,"l",nothing)
    Tc = params["Tc"]
    c1 = params["c1"]
    a  = epsilon_LorentzBerthelot(params["a"], k)
    b  = sigma_LorentzBerthelot(params["b"], l)

    epsilon_assoc = get!(params,"epsilon_assoc") do
        AssocParam("epsilon_assoc",_components)
    end

    bondvol = get!(params,"bondvol") do
        AssocParam("bondvol",_components)
    end

    bondvol,epsilon_assoc = assoc_mix(bondvol,epsilon_assoc,cbrt.(b),assoc_options,sites)
    packagedparams = CPAParam(Mw, Tc, a, b, c1, epsilon_assoc, bondvol)
    
    #init cubic model
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    init_alpha = init_model(alpha,components,alpha_userlocations,verbose)
    init_mixing = init_model(mixing,components,activity,mixing_userlocations,activity_userlocations,verbose)
    init_translation = init_model(translation,components,translation_userlocations,verbose)
    cubicparams = ABCubicParam(a, b, params["Tc"],Pc,Mw) #PR, RK, vdW
    init_cubicmodel = cubicmodel(_components,init_alpha,init_mixing,init_translation,cubicparams,init_idealmodel,String[])
    
    references = ["10.1021/ie051305v"]

    model = CPA(_components, radial_dist, init_cubicmodel, packagedparams, sites, init_idealmodel, assoc_options, references)
    set_reference_state!(model,reference_state;verbose)
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
    bondvol,epsilon_assoc = assoc_mix(bondvol,epsilon_assoc,cbrt.(b),assoc_options,model.sites) #combining rules for association

    model.params.epsilon_assoc.values.values[:] = epsilon_assoc.values.values
    model.params.bondvol.values.values[:] = bondvol.values.values
    return model
end

lb_volume(model::CPAModel,T,z) = lb_volume(model.cubicmodel,T,z)
T_scale(model::CPAModel,z) = T_scale(model.cubicmodel,z)

function p_scale(model::CPAModel,z)
    #does not depend on Pc, so it can be made optional on CPA input
    b = model.cubicmodel.params.b.values
    a = model.cubicmodel.params.a.values
    Ωa,Ωb = ab_consts(model.cubicmodel)
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
    return a_res(model.cubicmodel,V,T,z,_data) + a_assoc(model,V+c̄*n,T,z,_data)
end

ab_consts(model::CPAModel) = ab_consts(model.cubicmodel)

function Δ(model::CPAModel, V, T, z, i, j, a, b, _data = @f(data))
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
#optimized Δ function for CPA, we only calculate g once.
function  Δ(model::CPA, V, T, z,_data=@f(data))
    n,ā,b̄,c̄ = _data
    β = model.params.bondvol.values
    b_cubic = model.params.b.values
    η = n*b̄/(4*V)
    rdf = model.radial_dist
    g = if rdf == :CS #CPA original
        (1-0.5*η)/(1-η)^3
    elseif rdf == :KG #sCPA
        1/(1-1.9η)
    else
        zero(η)/zero(η)
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
