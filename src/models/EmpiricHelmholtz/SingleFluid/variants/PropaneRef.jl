#ancillary equations for calculation of P_sat, T_sat rhovsat y rholsat
#with the rework of single fluid models, now those functions are only used by other models, not PropaneRef.
function _propaneref_psat(T)
    T_c = 369.89
    P_c = 4.2512e6
    T>T_c && return zero(T)/zero(T)
    Tr = T/T_c
    θ = 1.0-Tr
    #    [-6.7722,1.6938,-1.3341,-3.1876,0.94937]
    #    [1.0,1.5,2.2,4.8,6.2]

    lnPsatPc = (-6.7722*θ + 1.6938*θ^1.5 -1.3341*θ^2.2 -3.1876*θ^4.8 + 0.94937*θ^6.2)/Tr
    Psat = exp(lnPsatPc)*P_c
    return Psat
end

function _propaneref_tsat(p)
    P_c = 4.2512e6
    T_c = 369.89
    p > P_c && return zero(p)/zero(p)
    #first aproximation
    A,B,C = 13.6515,1850.8,249.99-273.15
    T0 = B/(A - log(p*1e-3)) - C
    T0 > T_c && (T0 = T_c*p/P_c)
    f(T) = _propaneref_psat(T) - p
    prob = Roots.ZeroProblem(f,T0)
    return Roots.solve(prob,Roots.Order0())
end

"""
    PropaneRef <: EmpiricHelmholtzModel
    PropaneRef()

## Input parameters

None

## Description
Propane Reference Equation of State
```
δ = ρ/ρc
τ = T/Tc
a⁰(δ,τ) = log(δ) + n⁰₁ + n⁰₂τ + n⁰₃log(τ) + ∑n⁰ᵢ(1-exp(-γ⁰ᵢτ)), i ∈ 4:7
aʳ(δ,τ)  = aʳ₁+ aʳ₂ + aʳ₃
aʳ₁(δ,τ)  = ∑nᵢδ^(dᵢ)τ^(tᵢ), i ∈ 1:5
aʳ₂(δ,τ)  = ∑nᵢexp(-δ^cᵢ)δ^(dᵢ)τ^(tᵢ), i ∈ 6:11
aʳ₃(δ,τ)  = ∑nᵢexp(-ηᵢ(δ - εᵢ)^2 - βᵢ(τ - γᵢ)^2)δ^(dᵢ)τ^(tᵢ), i ∈ 12:18

```
Parameters  `n⁰`,`γ⁰`,`n`,`t`,`d`,`c`,`η`,`β`,`γ`,`ε` where obtained via fitting.

## References
1. Lemmon, E. W., McLinden, M. O., & Wagner, W. (2009). Thermodynamic properties of propane. III. A reference equation of state for temperatures from the melting line to 650 K and pressures up to 1000 MPa. Journal of Chemical and Engineering Data, 54(12), 3141–3180. [doi:10.1021/je900217v](https://doi.org/10.1021/je900217v)
"""
function PropaneRef()

    components = ["propane"]
    Mw = 44.09562 #g·mol-1
    T_c = 369.89    #K
    P_c = 4.2512e6 #Pa
    rho_c= 5000.0 # mol·m-3
    lb_volume = 1/53130
    Ttp = 85.525 #K
    ptp = 0.00017
    rhov_tp  = 2.4e-07
    rhol_tp = 16626.0
    Rgas = 8.314472
    acentric_factor = 0.1521

    properties = SingleFluidProperties(Mw,T_c,rho_c,lb_volume,T_c,P_c,rho_c,Ttp,ptp,rhov_tp,rhol_tp,acentric_factor,Rgas)

    a₁ = -4.970583
    a₂ = 4.29352
    u = -[1.062478, 3.344237,5.363757,11.762957]
    v = [3.043,5.874,9.337,7.922]
    c0 = 4.0 - 1

    ideal = SingleFluidIdealParam(a₁,a₂,c0,v,u)

    n = [0.042910051, 1.7313671, -2.4516524, 0.34157466, -0.46047898, -0.66847295, 0.20889705, 0.19421381, -0.22917851, -0.60405866, 0.066680654, 0.017534618, 0.33874242, 0.22228777, -0.23219062, -0.09220694, -0.47575718, -0.017486824]
    t = [1.0, 0.33, 0.8, 0.43, 0.9, 2.46, 2.09, 0.88, 1.09, 3.25, 4.62, 0.76, 2.5, 2.75, 3.05, 2.55, 8.4, 6.75]
    d = [4,1,1,2,2,1,3,6,6,2,3,1,1,1,2,2,4,1]
    l = [1,1,1,1,2,2]
    η = [0.963,1.977,1.917,2.307,2.546,3.28,14.6]
    β = [2.33,3.47,3.15,3.19,0.92,18.8,547.8]
    γ = [0.684,0.829,1.419,0.817,1.500,1.426,1.093]
    ε = [1.283,0.6936,0.788,0.473,0.8577,0.271,0.948]

    polexpgauss = PolExpGaussTerm(n,t,d,l,ones(length(l)),η,β,γ,ε)

    residual = SingleFluidResidualParam(;polexpgauss)

    anc_gas_fn = GenericAncEvaluator([-2.4887,-5.1069,-12.174,-30.495,-52.192,-134.89],[0.3785,1.07,2.7,5.5,10,20],T_c,rho_c,:exp,false)
    anc_liquid_fn = GenericAncEvaluator([1.82205,0.65802,0.21109,0.083973],[0.345,0.74,2.6,7.2],T_c,rho_c,:noexp,false)
    anc_ps_fn = GenericAncEvaluator([-6.7722,1.6938,-1.3341,-3.1876,0.94937],[1.0,1.5,2.2,4.8,6.2],T_c,P_c,:exp,true)
    ancillary_gas = PolExpVapour(anc_gas_fn)
    ancillary_liquid = PolExpLiquid(anc_liquid_fn)
    ancillary_pressure = PolExpSat(anc_ps_fn)
    ancillaries = CompositeModel(components,gas = ancillary_gas,liquid = ancillary_liquid,saturation = ancillary_pressure)

    references = ["1021/je900217v"]

    return SingleFluid(components,properties,ancillaries,ideal,residual,references)
end

function propane_ancillary_cs(components,T_c,P_c,Vc)
    rho_c = 1/Vc
    T_c0 = 369.89    #K
    P_c0 = 4.2512e6 #Pa
    rho_c0 = 5000.0 # mol·m-3

    anc_gas_fn = GenericAncEvaluator([-2.4887,-5.1069,-12.174,-30.495,-52.192,-134.89],[0.3785,1.07,2.7,5.5,10,20],T_c,rho_c,:exp,false)
    anc_liquid_fn = GenericAncEvaluator([1.82205,0.65802,0.21109,0.083973],[0.345,0.74,2.6,7.2],T_c,rho_c,:noexp,false)
    anc_ps_fn = GenericAncEvaluator([-6.7722,1.6938,-1.3341,-3.1876,0.94937],[1.0,1.5,2.2,4.8,6.2],T_c,P_c,:exp,true)
    ancillary_gas = PolExpVapour(anc_gas_fn)
    ancillary_liquid = PolExpLiquid(anc_liquid_fn)
    ancillary_pressure = PolExpSat(anc_ps_fn)
    ancillaries = CompositeModel(components,gas = ancillary_gas,liquid = ancillary_liquid,saturation = ancillary_pressure)
end

export PropaneRef