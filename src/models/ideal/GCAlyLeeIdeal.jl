struct GCAlyLeeParam <: EoSParam
    A::SingleParam{Float64}
    B::SingleParam{Float64}
    C::SingleParam{Float64}
    D::SingleParam{Float64}
    E::SingleParam{Float64}
    Mw::SingleParam{Float64}
    coeffs::SingleParam{NTuple{5,Float64}}
    reference_state::ReferenceState
end

abstract type GCAlyLeeModel <: IdealModel end
@newmodelgc GCAlyLeeIdeal GCAlyLeeModel GCAlyLeeParam false
default_references(::Type{GCAlyLeeIdeal}) = ["10.1021/acs.jced.5c00573"]
default_locations(::Type{GCAlyLeeIdeal}) = ["ideal/BurkhardtIdeal/GCAlyLeeIdeal.csv","properties/molarmass_groups.csv"]
default_gclocations(::Type{GCAlyLeeIdeal}) = ["ideal/BurkhardtIdeal/BurkhardtIdeal_groups.csv"]
default_ignore_missing_singleparams(::Type{GCAlyLeeIdeal}) = ["Mw"]

function transform_params(::Type{GCAlyLeeIdeal},params,groups)
    components = groups.components
    n = groups.n_flattenedgroups
    l = length(components)
    a,b,c,d,e = params["A"],params["B"],params["C"],params["D"],params["E"]
    _a,_b,_c,_d,_e = zeros(l),zeros(l),zeros(l),zeros(l),zeros(l)
    for i in 1:l
        #res +=z[i]*(log(z[i]/V))/Σz
        ni = n[i]
        _a[i] = (dot(a.values,ni) + 2.35963503)*8.31433/R̄
        _b[i] = (dot(b.values,ni) - 4.20519291)*8.31433/R̄
        _c[i] = (dot(c.values,ni) + 500.64232045)
        _d[i] = (dot(d.values,ni) + 3.44955031)*8.31433/R̄
        _e[i] = (dot(e.values,ni) + 514.21006282)
    end
    params["coeffs"] = reid_coeffs(_a,_b,_c,_d,_e,components)
    return params
end

function recombine_impl!(model::GCAlyLeeIdeal)
    coeffs = model.params.coeffs
    n = model.groups.n_flattenedgroups
    a = model.params.A.values
    b = model.params.B.values
    c = model.params.C.values
    d = model.params.D.values
    e = model.params.E.values
    for i in 1:length(model)
        #res +=z[i]*(log(z[i]/V))/Σz
        ni = n[i]
        _a = (dot(a,ni) + 2.35963503)*8.31433/R̄
        _b = (dot(b,ni) - 4.20519291)*8.31433/R̄
        _c = (dot(c,ni) + 500.64232045)
        _d = (dot(d,ni) + 3.44955031)*8.31433/R̄
        _e = (dot(e,ni) + 514.21006282)
        coeffs[i] = (_a,_b,_c,_d,_e)
    end
    return model
end

"""
    GCAlyLeeIdeal <: GCAlyLeeModel

    GCAlyLeeIdeal(components; 
    userlocations = String[],
    group_userlocations = String[]
    verbose = false)

## Input parameters

- `A`: Single Parameter (`Float64`)
- `B`: Single Parameter (`Float64`)
- `C`: Single Parameter (`Float64`)
- `D`: Single Parameter (`Float64`)
- `E`: Single Parameter (`Float64`)

## Description

AlyLee Group Contribution Ideal Model.

```
Cpᵢ(T)/R = A + B(CT⁻¹/sinh(CT⁻¹))² + D(ET⁻¹/cosh(ET⁻¹))²
A = ∑Nᵢₖ*(Aᵢₖ + 2.360)
B = ∑Nᵢₖ*(Bᵢₖ - 4.205)
C = ∑Nᵢₖ*(Cᵢₖ + 500.642)
D = ∑Nᵢₖ*(Dᵢₖ + 3.450)
E = ∑Nᵢₖ*(Eᵢₖ + 514.210)
```
## Group Fragmentation

Molecule fragmentation into functional groups is available in GCIdentifier.jl, using `Burkhardt2025Groups`

## References

1. Burkhardt, J., Bauer, G., Stierle, R., & Gross, J. (2026). A new group-contribution approach for ideal gas heat capacity, critical temperature and normal Boiling Point. Journal of Chemical and Engineering Data, 71(1), 6–23. [doi:10.1021/acs.jced.5c00573](http://doi.org/10.1021/acs.jced.5c00573)

"""
GCAlyLeeIdeal

export GCAlyLeeIdeal

function a_ideal(model::GCAlyLeeModel,V,T,z)
    #we transform from AlyLee terms to GERG2008 terms.
    coeffs = model.params.coeffs.values
    Σz = sum(z)
    res = zero(V+T+first(z))
    ρ = 1/V
    lnΣz = log(Σz)
    τi = one(T)/T
    logτi = log(τi)
    ng = model.groups.n_flattenedgroups
    @inbounds for i ∈ @comps
        #we suppose ρc = Tc = 1
        δi = ρ
        Ai,Bi,Ci,Di,Ei = coeffs[i]
        #integrate constant:
        #Tc,T0 = 1.0,298.15
        a₁,a₂,c₀ = _Cp0_constant_parse(Ai,1.0,298.15)
        #a₁ = Ai*-4.697596715569115 #(1 - log(τ0))
        #a₂ = -Ai*298.15 #T0/Tc
        c₀ -= 1
        ai = a₁ + a₂*τi + c₀*logτi
        zi = z[i]
        ai -= Di*EoSFunctions.logcosh(Ei*τi)
        ai += Bi*EoSFunctions.logabssinh(Ci*τi)
        res += xlogx(zi,δi)
        res += zi*ai
    end
    return res/Σz - 1
end