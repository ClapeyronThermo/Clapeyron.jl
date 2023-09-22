
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
    coeffs::SingleParam{NTuple{4,Float64}}
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

## List of available groups
|Name            |SMARTS                                              |Description    |
|----------------|----------------------------------------------------|---------------|
|-CH3           |[CX4H3]                                             |Methyl         |
|-CH2-           |[!R;CX4H2]                                          |Methylene      |
|>CH-            |[!R;CX4H]                                           |               |
|>C<             |[!R;CX4H0]                                          |               |
|CH2=CH-         |[CX3H2][CX3H1]                                      |               |
|-CH=CH-         |[CX3H1][CX3H1]                                      |               |
|=C<             |[\$([!R;#6X3H0]);!\$([!R;#6X3H0]=[#8])]               |               |
|=C=             |[\$([CX2H0](=*)=*)]                                  |               |
|CH              |[\$([CX2H1]#[!#7])]                                  |               |
|C               |[\$([CX2H0]#[!#7])]                                  |               |
|ring-CH2-       |[R;CX4H2]                                           |Cyclic alkane  |
|ring>CH-        |[R;CX4H]                                            |               |
|ring>C<         |[R;CX4H0]                                           |               |
|ring=CH-        |[R;CX3H1,cX3H1]                                     |Aromatic group |
|ring=C<         |[\$([R;#6X3H0]);!\$([R;#6X3H0]=[#8])]                 |               |
|-F             |[F]                                                 |Fluoride       |
|-Cl            |[Cl]                                                |Chloride       |
|-Br            |[Br]                                                |Bromide        |
|-I             |[I]                                                 |Iodide         |
|-OH (alcohol)  |[OX2H;!\$([OX2H]-[#6]=[O]);!\$([OX2H]-a)]             |Hydroxyl group |
|-OH (phenol)   |[O;H1;\$(O-!@c)]                                     |               |
|-O- (non-ring) |[OX2H0;!R;!\$([OX2H0]-[#6]=[#8])]                    |               |
|-O- (ring)     |[#8X2H0;R;!\$([#8X2H0]~[#6]=[#8])]                   |               |
|>C=O (non-ring) |[\$([CX3H0](=[OX1]));!\$([CX3](=[OX1])-[OX2]);!R]=O   |Ketone         |
|>C=O (ring)     |[\$([#6X3H0](=[OX1]));!\$([#6X3](=[#8X1])~[#8X2]);R]=O|               |
|O=CH- (aldehyde)|[CH;D2](=O)                                         |Aldehyde       |
|-COOH (acid)   |[OX2H]-[C]=O                                        |Carboxylic acid|
|-COO- (ester)  |[#6X3H0;!\$([#6X3H0](~O)(~O)(~O))](=[#8X1])[#8X2H0]  |Ester          |
|O (other than above)|[OX1H0;!\$([OX1H0]~[#6X3]);!\$([OX1H0]~[#7X3]~[#8])]  |Ketone         |
|-NH2           |[NX3H2]                                             |Amine          |
|>NH (non-ring)  |[NX3H1;!R]                                          |               |
|>NH (ring)      |[#7X3H1;R]                                          |               |
|>N- (non-ring)  |[#7X3H0;!\$([#7](~O)~O)]                             |               |
|-N= (non-ring) |[#7X2H0;!R]                                         |               |
|-N= (ring)     |[#7X2H0;R]                                          |               |
|=NH            |[#7X2H1]                                            |               |
|-CN            |[#6X2]#[#7X1H0]                                     |Nitrile        |
|-NO3           |[\$([#7X3,#7X3+][!#8])](=[O])~[O-]                   |Nitroxide      |
|-SH            |[SX2H]                                              |               |
|-S- (non-ring) |[#16X2H0;!R]                                        |               |
|-S- (ring)     |[#16X2H0;R]                                         |               |
"""
JobackIdeal

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

ReidIdeal(model::JobackIdeal) = ReidIdeal(model.components,ReidIdealParam(model.params.coeffs),model.references)

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

export JobackIdeal
