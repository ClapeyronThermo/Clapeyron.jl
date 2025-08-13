abstract type GrenkeElliottModel <: GibbsBasedModel end
@newmodelsingleton GrenkeElliottWater GrenkeElliottModel

"""
    GrenkeElliottWater <: GibbsBasedModel
    GrenkeElliottWater()

## Input parameters

None

## Description

Grenke and Elliott's model for liquid water at low temperatures (200-300 K) and high pressures (0.1−400 MPa)

```
v = MW*v₀*(1 - C*log((B + p)/(B + p₀)))
v₀ = a₁*exp(a₂T) + a₃*exp(a₄T) + a₅
p₀ = 101325
MW = 0.0180153
B = 1e8*(b₁/(1 + (T/b₂)^b₃)^b₄)
C = c₁/((1 + (T/c₂)^c₃)^c₄) + c₅
Cₚ(p₀,T) = (d₁exp(d₂T) + d₃)*Mw
g(p,T) = ∫∫Cₚ/T dT + ∫v dp + g₀₁ + g₀₂T
```

Where `g₀₁`and `g₀₂` are reference state constants, calculated to match `gibbs(ice,1 atm,273.15 K) == gibbs(liquid,1 atm,273.15 K)` and `Δh_fus(1 atm,273.15 K) = 6010 K J·mol⁻¹`

## References

1. Grenke, J. H., & Elliott, J. A. W. (2025). Analytic correlation for the thermodynamic properties of water at low temperatures (200-300 K) and high pressures (0.1-400 MPa). The Journal of Physical Chemistry. B, 129(7), 1997–2012. [doi:10.1021/acs.jpcb.4c03909](https://doi.org/10.1021/acs.jpcb.4c03909)
"""
GrenkeElliottWater

molecular_weight(model::GrenkeElliottModel,z) = 0.0180153*sum(z)

function v0_water(model::GrenkeElliottModel,T)
    a1 = 68.4089 #±24.0366 #[m3/kg]
    a2 = −0.0611145# ±0.0015535 #[1/K]
    a3 = 2.26928e-8# ±1.45185 × 10−8 #[m3/kg]
    a4 = 0.0215553# ±0.0018068 #[1/K]
    a5 = 9.88107e-4# ±0.01498 × 10−4 #[m3/kg]
    return a1*exp(a2*T) + a3*exp(a4*T) + a5
end

Base.length(model::GrenkeElliottModel) = 1

function B_water(model::GrenkeElliottModel,T)
    b1 = 3.1520397 #±1.0495712 [Pa]
    b2 = 203.8085375 #±96.2898756 [K]
    b3 = −11.1985548 #±4.6162330
    b4 = 7.2689427 #±33.5643396
    return 1e8*(b1/(1 + (T/b2)^b3)^b4)
end

function C_water(model::GrenkeElliottModel,T)
    c1 = 0.0790029 #±0.0949956
    c2 = 237.9009619 #±42.5843445 #[K]
    c3 = −14.8806681 #±12.8489102
    c4 = 0.8778057 #±2.7595783
    c5 = 0.0532605 #±0.1059223
    return c1/((1 + (T/c2)^c3)^c4) + c5
end

function volume_impl(model::GrenkeElliottModel,p,T,z,phase,threaded,vol0)
    #Tait-Taumann model
    mw = molecular_weight(model,z)
    v0 = v0_water(model,T)*mw
    B = B_water(model,T)
    C = C_water(model,T)
    p0 = 101325.0
    K1 = (B + p0)*exp(1/C)
    K2 = v0*C
    lnBp = log1p(p/B) - log1p(p0/B)
    return v0*(1 - C*lnBp)
end

function x0_pressure(model::GrenkeElliottModel,V,T,z)
    v = V/sum(z)
    mw = molecular_weight(model)
    v0 = v0_water(model,T)*mw
    B = B_water(model,T)
    C = C_water(model,T)
    p0 = 101325.0
    #=
    v = v0*(1 - C*lnBp)
    v/v0 = 1 - ClnBp
    ClnBp = 1 - v/v0
    =#
    p = B*expm1((1/C + log1p(p0/B)) - v/(v0*C))
    return p
end

function eos_g(model::GrenkeElliottModel,p,T,z)
    mw = molecular_weight(model,z)
    v0 = v0_water(model,T)*mw
    B = B_water(model,T)
    C = C_water(model,T)
    p0 = 101325.0
    K1 = (B + p0)*exp(1/C)
    K2 = v0*C
    lnBp = log1p(p/B) - log1p(p0/B)

    gT = water_g0(model,T)*mw #T-dependent
    gpTx = v0*(- C*(B + p)*lnBp + C*p + p) #p-T-dependent
    gpT0 = v0*(C*p0 + p0)
    gpT = gpTx - gpT0

    #=
    ## Reference state calculations ##

    We need to calculate g01 and g02 constants so that g(p,T) + g01 + g02*T is a valid water eos
    compatible with the ice equation of state.

    at T0 = 273.15, p0 = 101325:
    g(water) - g(ice) = 0       -> g_water - g_ice + g01 + g02*T0
    h(water) - h(ice) - dh = 0  -> h_water - h_ice - dh + h(g01 + g02*T0)

    h = g - T*∂g∂T
    h(g01 + g02*T) = g01 + g02*T - T*g02 = g01
    then:
    h(water) - h(ice) - c = 0
    h0 - h_ice - dh + g01 = 0
    g01 = dh + h_ice - h_water0

    g_water - g_ice + g01 + g02*T = 0
    g02 = (g_ice - g_water - g01)/T0

    Code to calculate constants:
    
    g_ice = 1.7703171209449624
    g_water0 = -94645.65737009811
    h_water0 = 20521.980322791598
    h_ice = -6005.5725370923665
    dh = 6010 #tabulated
    g01 = dh + h_ice - h_water0
    g02 = (g_ice - g_water0 - g01)/273.15
    =#
    g01, g02 = -20517.552859883966, 421.6180873040565
    gref = sum(z)*(g01 + g02*T)
    return gT + gpT + gref
end

function water_cp(model::GrenkeElliottModel,T)
    d1 = 4.44575e12 #± 3.00530e12 [J/kgK]
    d2 = −0.0928377 #± 0.0028242 [1/K]
    d3 = 4172.09 #± 8.06 [J/kgK]
    return d1*exp(d2*T) + d3
end

function water_g0(model::GrenkeElliottModel,T)
    #=
    for an expression of Cp = exp(k1*T), we solve:
        -T*d2gdT = exp(k*T)
    The solution is
        g(T) = -T*Ei(k*T) + exp(k*T)/k (+l1*T + l2)
    
    where Ei(x) is the exponential integral.
    =#
    d1 = 4.44575e12 #± 3.00530e12 [J/kgK]
    d2 = −0.0928377 #± 0.0028242 [1/K]
    d3 = 4172.09 #± 8.06 [J/kgK]

    #forgive me lord for this
    Ei = BlackBoxOptim.Distributions.SpecialFunctions.expinti

    exp_part = -T*Ei(d2*T) + exp(d2*T)/d2
    c_part = T - xlogx(T)
    return d1*exp_part + d3*c_part
end

export GrenkeElliottWater
