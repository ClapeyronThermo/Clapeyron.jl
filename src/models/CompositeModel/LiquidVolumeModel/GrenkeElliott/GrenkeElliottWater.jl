abstract type GrenkeElliottModel <: GibbsBasedModel end
@newmodelsingleton GrenkeElliottWater GrenkeElliottModel

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
    b1 = 3.340 #[Pa]
    b2 = 213.58 #[K]
    b3 = -11.252
    b4 = 4.490
    return 1e8*(b1/(1 + (T/b2)^b3)^b4)
end

function C_water(model::GrenkeElliottModel,T)
    c1 = 0.07169
    c2 = 225.18 # [K]
    c3 = -13.38
    c4 = 1.910
    c5 = 0.06829
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
    lnBp = log((B + p)/(B + p0))
    return v0*(1 - C*lnBp)
end

function eos_g(model::GrenkeElliottModel,p,T,z)
    mw = molecular_weight(model,z)
    v0 = v0_water(model,T)*mw
    B = B_water(model,T)
    C = C_water(model,T)
    p0 = 101325.0
    K1 = (B + p0)*exp(1/C)
    K2 = v0*C
    lnBp = log((B + p)/(B + p0))
    #V = v0*(1 - C*lnBp)
    gT = water_g0(model,T)*mw #T-dependent
    gpT = v0*(- C*(B + p)*lnBp + C*p + p) #p-T-dependent
    #=
    ## Reference state calculations ##

    We need to calculate l1 and l2 constants so that g(p,T) + l1 + l2*T is a valid water eos
    compatible with the ice equation of state.

    at triple point:
    g(water,Tt) - g(ice,Tt) = 0 -> = g_water - g_ice + l1 + l2*Tt
    h(water,Tt) - h(ice,Tt) - dh(Tt) = 0 h_water - h_ice - dh + h(l1 + l2*Tt)

    h = g - T*∂g∂T
    h(l1 + l2*T) = l1 + l2*T - T*l2 = l1
    then:
    h(water) - h(ice) - c = 0
    h0 - h_ice - dh + l1 = 0
    l1 = dh + h_ice - h_water0

    g_water - g_ice + l1 + l2*T = 0
    l2 = (g_ice - g_water - l1)/Triple
    g_ice = 0.011021455141630443
    g_water0 = -94649.62410803317
    h_water0 = 20522.78895968775
    h_ice = -6007.087598248805
    dh = 6007.098619700405 #calculated with IAPWS95 at triple point
    l1 = dh + h_ice - h_water0
    l2 = (g_ice - g_water0 - l1)/273.16
    =#
    l1, l2 = -20522.77793823615, 421.6298618674932
    gref = sum(z)*(l1 + l2*T)
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
    exp_part = -T*SpecialFunctions.expinti(d2*T) + exp(d2*T)/d2
    c_part = T - xlogx(T)
    return d1*exp_part + d3*c_part
end

#=
test function
function water_cp2(model,T)
    f(x) = water_g0(model,x)
    _,_,fx = Solvers.f∂f∂2f(f,T)
    return -T*fx
end
=#
