abstract type BerthelotModel <: vdWModel end

struct BerthelotParam <: EoSParam
    a::PairParam{Float64}
    b::PairParam{Float64}
    Tc::SingleParam{Float64}
    Pc::SingleParam{Float64}
    Vc::SingleParam{Float64}
    Mw::SingleParam{Float64}
end

struct Berthelot{T <: IdealModel,α,c,M} <: BerthelotModel
    components::Array{String,1}
    alpha::α
    mixing::M
    translation::c
    params::BerthelotParam
    idealmodel::T
    references::Array{String,1}
end

export Berthelot

"""
    Berthelot(components;
    idealmodel = BasicIdeal,
    alpha = ClausiusAlpha,
    mixing = vdW1fRule,
    activity = nothing,
    translation = NoTranslation,
    userlocations = String[],
    ideal_userlocations = String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    reference_state = nothing,
    verbose = false)

## Input parameters
- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Pressure `[Pa]`
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `Vc`: Single Parameter (`Float64`) - Molar Volume `[m^3/mol]`
- `k`: Pair Parameter (`Float64`) (optional)

## Model Parameters
- `a`: Pair Parameter (`Float64`)
- `b`: Pair Parameter (`Float64`)
- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Pressure `[Pa]`
- `Vc`: Single Parameter (`Float64`) - Molar Volume `[m^3/mol]`
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`


## Input models
- `idealmodel`: Ideal Model
- `alpha`: Alpha model
- `mixing`: Mixing model
- `activity`: Activity Model, used in the creation of the mixing model.
- `translation`: Translation Model

## Description

Berthelot Equation of state. it uses the Volume-Pressure Based mixing rules, that is:
```
a = 8*Pc*Vc^2
b = Vc/3
R = (8/3)*Pc*Vc/Tc
P = RT/(V-Nb) + a•α(T)/V²
α(T) = Tc/T
```

## References

1. Berthelot, D. (1899). Sur une méthode purement physique pour la détermination des poids moléculaires des gaz et des poids atomiques de leurs éléments. Journal de Physique Théorique et Appliquée, 8(1), 263–274. [doi:10.1051/jphystap:018990080026300](https://doi.org/10.1051/jphystap:018990080026300)

"""
Berthelot

function Berthelot(components;
    idealmodel = BasicIdeal,
    alpha = ClausiusAlpha,
    mixing = vdW1fRule,
    activity = nothing,
    translation = NoTranslation,
    userlocations = String[],
    ideal_userlocations = String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    verbose = false)
    formatted_components = format_components(components)
    params = getparams(components, ["properties/critical.csv", "properties/molarmass.csv","SAFT/PCSAFT/PCSAFT_unlike.csv"]; userlocations=userlocations, verbose=verbose)
    k  = get(params,"k",nothing)
    l = get(params,"l",nothing)
    pc = params["Pc"]
    Mw = params["Mw"]
    Tc = params["Tc"]
    Vc = params["Vc"]
    acentricfactor = get(params,"acentricfactor",nothing)
    init_mixing = init_model(mixing,components,activity,mixing_userlocations,activity_userlocations,verbose)
    a = PairParam("a",formatted_components,zeros(length(Tc)))
    b = PairParam("b",formatted_components,zeros(length(Tc)))
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose,reference_state)
    init_translation = init_model(translation,components,translation_userlocations,verbose)
    init_alpha = init_alphamodel(alpha,components,acentricfactor,alpha_userlocations,verbose)
    packagedparams = BerthelotParam(a,b,Tc,pc,Vc,Mw)
    references = String["10.1051/jphystap:018990080026300"]
    model = Berthelot(formatted_components,init_alpha,init_mixing,init_translation,packagedparams,init_idealmodel,references)
    recombine_cubic!(model,k,l)
end

function ab_premixing(model::BerthelotModel,mixing::MixingRule,k=nothing,l=nothing)
    #_Tc = model.params.Tc
    _pc = model.params.Pc
    _Vc = model.params.Vc
    a = model.params.a
    b = model.params.b
    diagvalues(a) .=  @. 3*_pc*_Vc^2
    diagvalues(b) .=  @. _Vc/3
    epsilon_LorentzBerthelot!(a,k)
    sigma_LorentzBerthelot!(b,l)
    #a = epsilon_LorentzBerthelot(SingleParam("a",components, @. Ωa*R̄^2*_Tc^2/_pc),kij)
    #b = sigma_LorentzBerthelot(SingleParam("b",components, @. Ωb*R̄*_Tc/_pc))
    return a,b
end

function pure_cubic_zc(model::BerthelotModel)
    return only(model.params.Pc.values)*only(model.params.Vc.values)/(R̄*only(model.params.Tc.values))
end

function a_res(model::BerthelotModel, V, T, z,_data = data(model,V,T,z))
    n,ā,b̄,c̄ = _data

    V̄c = 3*b̄
    #ā = 3*Pc*Vc^2*Tc
    P̄c = dot(model.params.Pc.values,z)/n
    T̄c = T*ā/(V̄c*V̄c*P̄c*3)
   # if T̄c isa Float64
   #     @show T̄c
   # end
    R = 8*P̄c*V̄c/(3*T̄c)
    #@show R
    #R = R̄
    RT⁻¹ = 1/(R*T)
    ρt = (V/n+c̄)^(-1) # translated density
    ρ  = n/V
    return (-log(1+(c̄-b̄)*ρ) - ā*ρt*RT⁻¹)
    #
    #return -log(V-n*b̄) - ā*n/(R̄*T*V) + log(V)
end

function T_scale(model::BerthelotModel,z = SA[1.0])
    comps = 1:length(z)
    Tc = model.params.Tc.values
    return sum(sqrt(Tc[i]*Tc[j])*z[i]*z[j] for i in comps for j in comps)/(sum(z)^2)
end

function p_scale(model::BerthelotModel,z = SA[1.0])
    dot(model.params.Pc.values,z)/sum(z)
end

function x0_crit_pure(model::BerthelotModel)
    lb_v = lb_volume(model)
    (1.1, log10(lb_v*3))
end

function crit_pure(model::BerthelotModel)
    single_component_check(crit_pure,model)
    Tc = model.params.Tc.values[1]
    Vc = model.params.Vc.values[1]
    Pc = pressure(model,Vc,Tc)
    #return Base.invoke(crit_pure,Tuple{EoSModel},model)
    return (Tc,Pc,Vc)
end


#=
z = PV/RT = 1/(1 - b/v) - a/RT2V2
P = RT/(v - b)  - a/TV

=#