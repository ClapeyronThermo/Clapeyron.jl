
#TODO
# - standardize chemical group notation to make it the same as SAFTgammaMie
# - add a database of group Mw

struct JobackIdealParam <: EoSParam
    Mw_gc::SingleParam{Float64}
    Mw::SingleParam{Float64}
    N_a::SingleParam{Int}
    T_c::SingleParam{Float64}
    P_c::SingleParam{Float64}
    V_c::SingleParam{Float64}
    T_b::SingleParam{Float64}
    T_m::SingleParam{Float64}
    H_form::SingleParam{Float64}
    G_form::SingleParam{Float64}
    a::SingleParam{Float64}
    b::SingleParam{Float64}
    c::SingleParam{Float64}
    d::SingleParam{Float64}
    H_fusion::SingleParam{Float64}
    H_vap::SingleParam{Float64}
    eta_a::SingleParam{Float64}
    eta_b::SingleParam{Float64}
    coeffs::SingleParam{NTuple{5,Float64}}
    reference_state::ReferenceState
end

abstract type JobackIdealModel <: ReidIdealModel end
@newmodelgc JobackIdeal JobackIdealModel JobackIdealParam false
default_references(::Type{JobackIdeal}) = ["10.1080/00986448708960487"]
default_locations(::Type{JobackIdeal}) = ["ideal/JobackIdeal.csv","properties/molarmass_groups.csv","properties/natoms_groups.csv"]
default_gclocations(::Type{JobackIdeal}) = ["ideal/JobackIdeal_Groups.csv"]
function transform_params(::Type{JobackIdeal},params,groups)
    components = groups.components
    n = groups.n_flattenedgroups
    l = length(components)
    i_groups = groups.i_groups
    a,b,c,d = params["a"],params["b"],params["c"],params["d"]
    Mw_gc = params["Mw"]
    params["Mw_gc"] = Mw_gc
    params["Mw"] = SingleParam("Mw",components,group_Mw(Mw_gc,groups))
    _a,_b,_c,_d = zeros(l),zeros(l),zeros(l),zeros(l)
    for i in 1:l
        #res +=z[i]*(log(z[i]/V))/Σz
        ni = n[i]
        _a[i] = ∑(a[j]*ni[j] for j in i_groups[i]) - 37.93
        _b[i] = ∑(b[j]*ni[j] for j in i_groups[i]) + 0.210
        _c[i] = ∑(c[j]*ni[j] for j in i_groups[i]) - 3.91e-4
        _d[i] = ∑(d[j]*ni[j] for j in i_groups[i]) + 2.06e-7
    end
    params["coeffs"] = reid_coeffs(_a,_b,_c,_d,components)
    return params
end

export JobackIdeal

"""
    JobackIdeal <: JobackIdealModel

    JobackIdeal(components; 
    userlocations = String[],
    reference_state = nothing,
    verbose = false)

## Input parameters

- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `N_a`: Single Parameter (`Float64`)
- `T_c`: Single Parameter (`Float64`)
- `P_c`: Single Parameter (`Float64`)
- `V_c`: Single Parameter (`Float64`)
- `T_b`: Single Parameter (`Float64`)
- `T_m`: Single Parameter (`Float64`)
- `H_form`: Single Parameter (`Float64`)
- `G_form`: Single Parameter (`Float64`)
- `a`: Single Parameter (`Float64`)
- `b`: Single Parameter (`Float64`)
- `c`: Single Parameter (`Float64`)
- `d`: Single Parameter (`Float64`)
- `H_fusion`: Single Parameter (`Float64`)
- `H_vap`: Single Parameter (`Float64`)
- `eta_a`: Single Parameter (`Float64`)
- `eta_b`: Single Parameter (`Float64`)

## Description
 
Joback Group Contribution Ideal Model. GC version of `ReidIdeal`. Helmholtz energy obtained via integration of specific heat capacity:

```
aᵢ = ∑(νᵢₖbₖ) - 37.93
bᵢ = ∑(νᵢₖbₖ) + 0.210
cᵢ = ∑(νᵢₖcₖ) - 3.91e-4
dᵢ = ∑(νᵢₖbₖ) + 2.06e-7
Cpᵢ(T) = aᵢ  + bᵢT + cᵢT^2 + dᵢT^3
```

The GC-averaged Reid Model is available by using `ReidIdeal(model::JobackIdeal)`.

The estimated critical point of a single component can be obtained via `crit_pure(model::JobackIdeal)`
## References
1. Joback, K. G., & Reid, R. C. (1987). Estimation of pure-component properties from group-contributions. Chemical Engineering Communications, 57(1–6), 233–243. [doi:10.1080/00986448708960487](https://doi.org/10.1080/00986448708960487)

## List of available groups
|Name    |Description         |
|--------|--------------------|
|-CH3   |Methyl              |
|-CH2-   |Methylene           |
|>CH-    |                    |
|>C<     |                    |
|CH2=CH- |                    |
|-CH=CH- |                    |
|=C<     |                    |
|=C=     |                    |
|CH      |                    |
|C       |                    |
|ring-CH2-|Cyclic alkane       |
|ring>CH-|                    |
|ring>C< |                    |
|ring=CH-|Aromatic group      |
|ring=C< |                    |
|-F     |Fluoride            |
|-Cl    |Chloride            |
|-Br    |Bromide             |
|-I     |Iodide              |
|-OH (alcohol)|Hydroxyl group      |
|-OH (phenol)|                    |
|-O- (non-ring)|                    |
|-O- (ring)|                    |
|>C=O (non-ring)|Ketone              |
|>C=O (ring)|                    |
|O=CH- (aldehyde)|Aldehyde            |
|-COOH (acid)|Carboxylic acid     |
|-COO- (ester)|Ester               |
|O (other than above)|Ketone              |
|-NH2   |Amine               |
|>NH (non-ring)|                    |
|>NH (ring)|                    |
|>N- (non-ring)|                    |
|-N= (non-ring)|                    |
|-N= (ring)|                    |
|=NH    |                    |
|-CN    |Nitrile             |
|-NO3   |Nitroxide           |
|-SH    |                    |
|-S- (non-ring)|                    |
|-S- (ring)|                    |
"""
JobackIdeal
mw(model::JobackIdeal) = model.params.Mw.values
molecular_weight(model::JobackIdeal,z) = molecular_weight(ReidIdeal(model),z)
function recombine_impl!(model::JobackIdeal)
    coeffs = model.params.coeffs
    i_groups = model.groups.i_groups
    n = model.groups.n_flattenedgroups
    a = model.params.a.values
    b = model.params.b.values
    c = model.params.c.values
    d = model.params.d.values
    for i in 1:length(model)
        #res +=z[i]*(log(z[i]/V))/Σz
        ni = n[i]
        _a = ∑(a[j]*ni[j] for j in i_groups[i]) - 37.93
        _b = ∑(b[j]*ni[j] for j in i_groups[i]) + 0.210
        _c = ∑(c[j]*ni[j] for j in i_groups[i]) - 3.91e-4
        _d = ∑(d[j]*ni[j] for j in i_groups[i]) + 2.06e-7
        coeffs[i] = (_a,_b,_c,_d)
    end
    return model
end

function ReidIdeal(model::JobackIdeal)
    comps = model.components
    coeffs = model.params.coeffs
    Mw = model.params.Mw,model
    a = SingleParam("a",comps,getindex.(coeffs.values,1))
    b = SingleParam("b",comps,getindex.(coeffs.values,2))
    c = SingleParam("c",comps,getindex.(coeffs.values,3))
    d = SingleParam("d",comps,getindex.(coeffs.values,4))
    e = SingleParam("e",comps,getindex.(coeffs.values,5))
    reference_state = model.params.reference_state
    param = ReidIdealParam(a,b,c,d,e,coeffs,reference_state,Mw)
    ReidIdeal(model.components,param,model.references)
end

"""
    JobackGC

Module containing group contribution calculations using the joback method. the available functions are:

- `JobackGC.T_c(model::JobackModel)`: critical temperature (in K)
- `JobackGC.P_c(model::JobackModel)`: critical pressure (in Pa)
- `JobackGC.V_c(model::JobackModel)`: critical volume (in m3/mol)
- `JobackGC.T_b(model::JobackModel)`: normal boiling point (in K)
- `JobackGC.H_form(model::JobackModel)`: enthalpy of formation at 298K, ideal gas (in J/mol)
- `JobackGC.G_form(model::JobackModel)`: gibbs energy of formation at 298K, ideal gas (in J/mol)
- `JobackGC.S_form(model::JobackModel)`: entropy of formation at 298K, ideal gas (in J/mol/K)
- `JobackGC.H_fusion(model::JobackModel)`: enthalpy of fusion (in J/mol, at 1 atm)
- `JobackGC.H_vap(model::JobackModel)`: molar enthalpy of vaporization (in J/mol, at normal boiling point)
- `JobackGC.C_p(model::JobackModel, T)`: ideal gas isobaric heat capacity (in J/mol/K)
- `JobackGC.Visc(model::JobackModel, T)`: liquid dynamic viscocity (in Pa*s)
"""
module JobackGC
    using Clapeyron: JobackIdeal
    using LinearAlgebra: dot
"""
    JobackGC.T_b(model::JobackIdeal)::Vector{Float64}

Given a `JobackIdeal` model, returns a vector containing normal boiling points (in K), for each component.
"""
function T_b(model::JobackIdeal)
    result = zeros(Float64,length(model))
    n = model.groups.n_flattenedgroups
    T_b = model.params.T_b.values
    for i in 1:length(model)
        result[i] = dot(n[i],T_b) + 198.2
    end
    return result
end

"""
    JobackGC.T_c(model::JobackIdeal)::Vector{Float64}

Given a `JobackIdeal` model, returns a vector containing critical temperatures (in K), for each component.
"""
function T_c(model::JobackIdeal,Tb=T_b(model))
    n = model.groups.n_flattenedgroups
    T_c = model.params.T_c.values
    result = zeros(Float64,length(model))
    for i in 1:length(model)
        ΣT_ci = dot(n[i],T_c)
        result[i] = Tb[i] / (0.584+0.965*ΣT_ci - ΣT_ci^2)
    end
    return result
end

"""
    JobackGC.V_c(model::JobackIdeal)::Vector{Float64}

Given a `JobackIdeal` model, returns a vector containing critical volumes (in m3/mol), for each component.
"""
function V_c(model::JobackIdeal)
    n = model.groups.n_flattenedgroups
    V_c = model.params.V_c.values
    result = zeros(Float64,length(model))
    for i in 1:length(model)
        result[i] = dot(n[i],V_c) + 17.5
    end
    result .*= 1e-6
    return result
end

"""
    JobackGC.P_c(model::JobackIdeal)::Vector{Float64}

Given a `JobackIdeal` model, returns a vector containing critical pressures (in Pa), for each component.
"""
function P_c(model::JobackIdeal)
    n = model.groups.n_flattenedgroups
    P_c = model.params.P_c.values
    N_a = model.params.N_a.values
    result = zeros(Float64,length(model))
    for i in 1:length(model)
        ΣP_ci  = dot(n[i],P_c)
        ΣN_a = dot(n[i],N_a)
        result[i] = (0.113 + 0.0032*ΣN_a - ΣP_ci)^-2
    end
    #result in bar, converting to Pa
    result .*= 1e5
    return result
end

"""
    G_form(model::JobackIdeal)::Vector{Float64}

Given a `JobackIdeal` model, returns a vector containing gibbs energies of formation (J/mol, at ideal gas, 298K), for each component.
"""
function G_form(model::JobackIdeal)
    n = model.groups.n_flattenedgroups
    G_form = model.params.G_form.values
    result = zeros(Float64,length(model))
    for i in 1:length(model)
        result[i] = dot(n[i],G_form) + 53.88
    end
    #result in kJ, converting to J
    result .*= 1000
    return result
end

"""
    H_form(model::JobackIdeal)::Vector{Float64}

Given a `JobackIdeal` model, returns a vector containing entalpies of formation (J/mol, at ideal gas, 298K), for each component.
"""
function H_form(model::JobackIdeal)
    n = model.groups.n_flattenedgroups
    H_form = model.params.H_form.values
    result = zeros(Float64,length(model))
    for i in 1:length(model)
        result[i] = dot(n[i],H_form) + 68.29
    end
    #result in kJ, converting to J
    result .*= 1000
    return result
end

"""
    S_form(model::JobackIdeal)::Vector{Float64}

Given a `JobackIdeal` model, returns a vector containing entropies of formation (J/mol/K, at ideal gas, 298K), for each component.
"""
function S_form(model::JobackIdeal)
    n = model.groups.n_flattenedgroups
    H_form = model.params.H_form.values
    G_form = model.params.G_form.values
    result = zeros(Float64,length(model))
    for i in 1:length(model)
        H_formᵢ = dot(n[i],H_form) + 68.29
        G_formᵢ = dot(n[i],G_form) + 53.88
        result[i] = (H_formᵢ - G_formᵢ)/298
    end
    #result in kJ, converting to J
    result .*= 1000
    return result
end

"""
    H_fusion(model::JobackIdeal)::Vector{Float64}

Given a `JobackIdeal` model, returns a vector containing enthalpies of fusion (J/mol, at 1 atm), for each component.
"""
function H_fusion(model::JobackIdeal)
    n = model.groups.n_flattenedgroups
    H_fusion = model.params.H_fusion.values
    result = zeros(Float64,length(model))
    for i in 1:length(model)
        result[i] = dot(n[i],H_fusion) - 0.88
    end
    #result in kJ, converting to J
    result .*= 1000
    return result
end

"""
    H_vap(model::JobackIdeal)::Vector{Float64}

Given a `JobackIdeal` model, returns a vector containing enthalpies of vaporization (J/mol, at normal boiling point), for each component.
"""
function H_vap(model::JobackIdeal)
    n = model.groups.n_flattenedgroups
    H_vap = model.params.H_vap.values
    result = zeros(Float64,length(model))
    for i in 1:length(model)
        result[i] = dot(n[i],H_vap) + 15.3
    end
    #result in kJ, converting to J
    result .*= 1000
    return result
end

"""
    C_p(model::JobackIdeal,T)::Vector

Given a `JobackIdeal` model, returns a vector containing ideal gas isobaric heat capacities (J/mol/K), for each component.
"""
function C_p(model::JobackIdeal,T)
    n = model.groups.n_flattenedgroups
    a = model.params.a.values
    b = model.params.b.values
    c = model.params.c.values
    d = model.params.d.values
    result = zeros(Base.promote_eltype(T,Float64),length(model))
    for i in 1:length(model)
        ai = dot(n[i],a) - 37.93
        bi = dot(n[i],b) + 0.210
        ci = dot(n[i],c) - 3.91e-4
        di = dot(n[i],d) + 2.06e-7
        result[i] = evalpoly(T, (ai,bi,ci,di))
    end
    return result
end

"""
    Visc(model::JobackIdeal,T)::Vector

Given a `JobackIdeal` model, returns a vector containing liquid dynamic viscocities (Pa*s), for each component.
"""
function Visc(model::JobackIdeal,T)
    n = model.groups.n_flattenedgroups
    ηa = model.params.eta_a.values
    ηb = model.params.eta_b.values
    Mw = model.params.Mw_gc.values
    result = zeros(Base.promote_eltype(T,Float64),length(model))
    for i in 1:length(model)
        ηai = dot(n[i],ηa)
        ηbi = dot(n[i],ηb)
        Mwi = dot(n[i],Mw)
        result[i] = Mwi*exp((ηai - 597.82)/T + ηbi - 11.202)
    end
    return result
end

end #JobackGC module

function crit_pure(model::JobackIdeal)
    return (JobackGC.T_c(model)[1],JobackGC.P_c(model)[1],JobackGC.V_c(model)[1])
end

export JobackIdeal, JobackGC
