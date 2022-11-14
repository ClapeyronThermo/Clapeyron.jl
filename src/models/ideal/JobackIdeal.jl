
#TODO
# - standarize chemical group notation to make it the same as SAFTgammaMie
# - add a database of group Mw

struct JobackIdealParam <: EoSParam
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
end


abstract type JobackIdealModel <: IdealModel end

struct JobackIdeal <: JobackIdealModel
    components::Array{String,1}
    groups::GroupParam
    params::JobackIdealParam
    reidmodel::ReidIdeal
    references::Array{String,1}
end
@registermodel JobackIdeal

export JobackIdeal

"""
    JobackIdeal <: JobackIdealModel
    JobackIdeal(components; 
    userlocations::Array{String,1}=String[], 
    verbose=false)

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

"""
JobackIdeal

function JobackIdeal(components;userlocations=String[], verbose=false, kwargs...)
    groups = GroupParam(components,["ideal/JobackIdeal_Groups.csv"], verbose=verbose)
    params = getparams(groups, ["ideal/JobackIdeal.csv","properties/molarmass_groups.csv"]; userlocations=userlocations, verbose=verbose)
    Mw = params["Mw"]::SingleParam{Float64}
    N_a = params["N_a"]::SingleParam{Int}
    T_c = params["T_c"]::SingleParam{Float64}
    P_c = params["P_c"]::SingleParam{Float64}
    V_c = params["V_c"]::SingleParam{Float64}
    T_b = params["T_b"]::SingleParam{Float64}
    T_m = params["T_m"]::SingleParam{Float64}
    H_form = params["H_form"]::SingleParam{Float64}
    G_form = params["G_form"]::SingleParam{Float64}
    a = params["a"]::SingleParam{Float64}
    b = params["b"]::SingleParam{Float64}
    c = params["c"]::SingleParam{Float64}
    d = params["d"]::SingleParam{Float64}
    H_fusion = params["H_fusion"]::SingleParam{Float64}
    H_vap = params["H_vap"]::SingleParam{Float64}
    eta_a = params["eta_a"]::SingleParam{Float64}
    eta_b = params["eta_b"]::SingleParam{Float64}
    packagedparams = JobackIdealParam(
    Mw,
    N_a,
    T_c,
    P_c,
    V_c,
    T_b,
    T_m,
    H_form,
    G_form,
    a,
    b,
    c,
    d,
    H_fusion,
    H_vap,
    eta_a,
    eta_b)
    references = ["10.1080/00986448708960487"]
    
    comps = 1:length(groups.components)
    i_groups = groups.i_groups
    n = groups.n_flattenedgroups
    coeffs = Vector{NTuple{4,Float64}}(undef,length(comps)) 
    reidparam = ReidIdealParam(SingleParam("GC-averaged Reid Coefficients",groups.components,coeffs))
    reidmodel = ReidIdeal(reidparam)
    model = JobackIdeal(groups.components,groups,packagedparams,reidmodel,references)
    recombine!(model)
    return model
end

function recombine_impl!(model::JobackIdeal)
    coeffs = model.reidmodel.params.coeffs
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

ReidIdeal(model::JobackIdeal) = model.reidmodel

function VT_isobaric_heat_capacity(model::JobackIdeal,V,T,z=SA[1.])
    return VT_isobaric_heat_capacity(model.reidmodel,V,T,z)
end

function T_b(model::JobackIdeal)
    n = model.groups.n_flattenedgroups
    T_b = model.params.T_b.values
    res = 0.0
    @inbounds begin
        ni = n[1]
        res += ∑(T_b[j]*ni[j] for j in @groups(1))
    end
    res = res + 198.2
end

function T_c(model::JobackIdeal,Tb=T_b(model))
    n = model.groups.n_flattenedgroups
    T_c = model.params.T_c.values
    res = 0.0
    @inbounds begin
        ni = n[1]
        ΣT_ci = ∑(T_c[j]*ni[j] for j in @groups(1))
        res += Tb * (0.584+0.965*ΣT_ci - ΣT_ci^2)^-1
    end
    return res
end

#this style of writing is uglier but faster,
#ideally they should be of the same speed
function V_c(model::JobackIdeal)
    n = model.groups.n_flattenedgroups
    V_c = model.params.V_c.values
    res = 0.0
    groups = @groups(1)
    ni = n[1]
    @inbounds begin
        ΣV_ci = zero(res)
        for idx in 1:length(groups)
            j = groups[idx]
            ΣV_ci +=V_c[j]*ni[j]
        end
        res += 17.5 + ΣV_ci
    end
    #result in mL/mol, converting to m3/mol
    return res*1e-6
end

function P_c(model::JobackIdeal)
    n = model.groups.n_flattenedgroups
    P_c = model.params.P_c.values
    N_a = model.params.N_a.values
    groups = @groups(1)
    ni = n[1]
    @inbounds begin
        ΣP_ci  = 0.0
        ΣN_a = 0
        for idx in 1:length(groups)
            j = groups[idx]
            ΣP_ci += P_c[j]*ni[j]
            ΣN_a += N_a[j]*ni[j]
        end
        res = (0.113 + 0.0032*ΣN_a - ΣP_ci)^-2
    end
    #result in bar, converting to Pa
    return res*100000.0
end

function crit_pure(model::JobackIdeal)
    return (T_c(model),P_c(model),V_c(model))
end

function a_ideal(model::JobackIdealModel, V, T, z)
    return a_ideal(model.reidmodel,V,T,z)
end

export JobackIdeal


##utilities
#=
function reid_to_joback(a,b,c,d)
    _a = a + 37.93
    _b = b - 0.210
    _c = c + 3.91e-4
    _d = d - 2.06e-7
    return (_a,_b,_c,_d)
end

#Vc = 5.594807453383915e-5
#Tc = 647.096
#Pc = 2.2064e7
#Tb = 373.15
function crit_to_joback(Tb,Tc,Pc,Vc,na)
    tb = Tb - 198.2
    vc = Vc*1e6 - 17.5
    Tx = Tc/Tb
    Tx = 1/Tx - 0.584 #(0.965Ti(1 - Ti))
    Tx = Tx/0.965
    #Tx = Ti(1-Ti), #-Ti2 + Ti - Tx = 0, Ti2 - Ti + Tx = 0 
    tc = (1 + sqrt(1 -4Tx))/2
    fxt(t) = Tb*(0.584 + 0.965*t - t*t)^-1 - Tc
    tc = Roots.find_zero(fxt,tc)
    fxp(p) = (0.113 + 0.0032*na - p)^-2 - 1e-5*Pc
    pc = sqrt(abs(1/(Pc*Pc*1e-10) -0.113 - 0.0032*na)) 
    @show pc
    pc = Roots.find_zero(fxp,pc)
    return (tc,pc,vc,tb)
end
=#
