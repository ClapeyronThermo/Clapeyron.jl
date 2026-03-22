struct BurkhardtIdealParam <: EoSParam
    A::SingleParam{Float64}
    B::SingleParam{Float64}
    C::SingleParam{Float64}
    D::SingleParam{Float64}
    Mw::SingleParam{Float64}
    reference_state::ReferenceState
end

abstract type BurkhardtIdealModel <: IdealModel end
@newmodelgc BurkhardtIdeal BurkhardtIdealModel BurkhardtIdealParam false
default_references(::Type{BurkhardtIdeal}) = ["10.1021/acs.jced.5c00573"]
default_locations(::Type{BurkhardtIdeal}) = ["ideal/BurkhardtIdeal/BurkhardtIdeal.csv","properties/molarmass_groups.csv"]
default_gclocations(::Type{BurkhardtIdeal}) = ["ideal/BurkhardtIdeal/BurkhardtIdeal_groups.csv"]
default_ignore_missing_singleparams(::Type{BurkhardtIdeal}) = ["Mw"]

"""
    BurkhardtIdeal <: BurkhardtIdealModel

    BurkhardtIdeal(components; 
    userlocations = String[],
    group_userlocations = String[]
    verbose = false)

## Input parameters

- `A`: Single Parameter (`Float64`)
- `B`: Single Parameter (`Float64`)
- `C`: Single Parameter (`Float64`)
- `D`: Single Parameter (`Float64`)

## Description

Burkhardt Group Contribution Ideal Model. Based on the AlyLee model, with a different summation rule
```
Cpᵢ(T)/R = ∑νᵢₖAₖ + c₁ + ∑[Bᵢₖ(CᵢₖT⁻¹/sinh(CᵢₖT⁻¹))² + Dᵢₖ(EᵢₖT⁻¹/cosh(EᵢₖT⁻¹))²]
Bᵢₖ = Nᵢₖ*(Bₖ + c₂)
Cᵢₖ = Cₖ + c₃
Dᵢₖ = Nᵢₖ*(Dₖ + c₄)
Eᵢₖ = c₅*(Cₖ + c₃)
c₁,c₂,c₃,c₄,c₅ = -0.350,2.245,8.858,1.487,0.432
```

## Group Fragmentation

Molecule fragmentation into functional groups is available in GCIdentifier.jl, using `Burkhardt2025Groups`

## References

1. Burkhardt, J., Bauer, G., Stierle, R., & Gross, J. (2026). A new group-contribution approach for ideal gas heat capacity, critical temperature and normal Boiling Point. Journal of Chemical and Engineering Data, 71(1), 6–23. [doi:10.1021/acs.jced.5c00573](http://doi.org/10.1021/acs.jced.5c00573)

"""
BurkhardtIdeal

export BurkhardtIdeal

function a_ideal(model::BurkhardtIdealModel,V,T,z)
    #we transform from AlyLee terms to GERG2008 terms.

    A = model.params.A.values
    B = model.params.B.values
    C = model.params.C.values
    D = model.params.D.values
    
    Σz = sum(z)
    res = zero(V+T+first(z))
    ρ = 1/V
    lnΣz = log(Σz)
    τi = one(T)/T
    logτi = log(τi)
    ng = model.groups.n_flattenedgroups
    c1,c2,c3,c4,c5 = -0.350,2.245,8.858,1.487,0.432
    @inbounds for i ∈ @comps
        #we suppose ρc = Tc = 1
        δi = ρ
        Ni = ng[i]
        Ai = (dot(Ni,A) + c1)*8.31433/R̄
        #integrate constant:
        #Tc,T0 = 1.0,298.15
        a₁,a₂,c₀ = _Cp0_constant_parse(Ai,1.0,298.15)
        #a₁ = Ai*-4.697596715569115 #(1 - log(τ0))
        #a₂ = -Ai*298.15 #T0/Tc
        c₀ -= 1
        ai = a₁ + a₂*τi + c₀*logτi
        zi = z[i]

        for α in model.groups.i_groups[i]
            Niα = Ni[α]
            Cα = C[α]
            Biα = Niα*(B[α] + c2)*8.31433/R̄
            Ciα = 100*(Cα + c3)
            Diα = Niα*(D[α] + c4)*8.31433/R̄
            Eiα = Ciα*c5
            ai -= Diα*EoSFunctions.logcosh(Eiα*τi)
            ai += Biα*EoSFunctions.logabssinh(Ciα*τi)
        end

        res += xlogx(zi,δi)
        res += zi*ai
    end
    return res/Σz - 1
end
