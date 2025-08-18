abstract type HoltenWaterModel <: GibbsBasedModel end
@newmodelsingleton HoltenWater HoltenWaterModel

"""
    GrenkeElliottWater <: GibbsBasedModel
    GrenkeElliottWater()

## Input parameters

None

## Description

Holten's model for liquid water at low temperatures (homogeneous ice nucleation temperature up to
300 K) and high pressures (up to 400 MPa, but can be extrapolated up to 1000 MPa).

```
g(P,T) = R*Tc*ĝ
ĝ = ĝᴬ + T̂*(x*L + x*log(x) + (1 - x)*log(1 - x) + ω*x*(1-x))
x: solution of L + log(x) - log(1 - x) + ω*(1 - 2*x) = 0
L = L₀*K₂*(1 + k₀*k₂ + k₁*(p + k₂*t) - K₁)/(2*k₁*k₂)
K₁ = sqrt((1 + k₀*k₂ + k₁*(p - k₂*t))^2 - 4*k₀*k₁*k₂*(p - k₂*t))
K₂ = sqrt(1 + k₂*k₂)
t = (T - Tc)/Tc
p = P/(R*Tc*ρ0)
ω = 2 + ω₀*p
ĝᴬ = ∑(cᵢ * τ^aᵢ * π^bᵢ * exp(-dᵢ*π))
```

Where x is the fraction of water molecules in with low-density structure.

## References

1. Holten, V., Sengers, J. V., & Anisimov, M. A. (2014). Equation of state for supercooled water at pressures up to 400 MPa. Journal of Physical and Chemical Reference Data, 43(4), 043101. [doi:10.1063/1.4895593](https://doi.org/10.1063/1.4895593)
"""
HoltenWater

function eos_g(model::HoltenWaterModel,p,T,z)
    Tc = 228.2
    Rm = 461.523087
    ρ0 = 1081.6482
    𝕡 = p/(Rm*Tc*ρ0)
    L = water_L(model,𝕡,T)
    ω = 2 + 0.5212269*𝕡
    x = water_x_frac(model,L,ω)
    
    T̂ = T/Tc
    ĝᴬ = water_eos_g_a(model,p,T)
    ĝ = ĝᴬ + T̂*(x*L + xlogx(x) + xlogx(1-x) + ω*x*(1-x))
    g = ĝ*Rm*Tc*molecular_weight(model,z)
    return g
end

molecular_weight(model::HoltenWaterModel,z) = 0.0180153*sum(z)
default_references(::Type{HoltenWater}) = ["10.1063/1.4895593"]

function water_L(model::HoltenWaterModel,p,T)
    L0 = 0.76317954
    k0 = 0.072158686
    k₁ = -0.31569232
    k₂ = 5.2992608
    Tc = 228.2
    t = (T - Tc)/Tc
    K₁ = sqrt((1 + k0*k₂ + k₁*(p - k₂*t))^2 - 4*k0*k₁*k₂*(p - k₂*t))
    K₂ = sqrt(1 + k₂*k₂)
    L = L0*K₂*(1 + k0*k₂ + k₁*(p + k₂*t) - K₁)/(2*k₁*k₂)
end

function water_x_frac(model::HoltenWaterModel,L,ω)
    
    k0 = exp(-L - ω)
    function f(x)
            
        #fx = L + log(x) - log(1 - x) + ω*(1 - 2*x)
        #dfx = 1/x + 1/(1 - x) - 2*x*ω
        #d2fx = -1/(x*x) - 1/((1 - x)*(1 - x)) - 2*ω
        #@show fx
        #return fx,fx/dfx,d2fx/dfx
        k = k0*exp(2*ω*x)
        xx = k/(1 + k)
        fx = xx - x
        dfx = 2*ω*xx/(1 + k) - 1
        return fx,fx/dfx
    end
    prob = Roots.ZeroProblem(f,k0)
    return Roots.solve(prob,Roots.Newton())
end

function water_eos_g_a(model::HoltenWaterModel,p,T)
    c = HoltenWaterConsts.c
    a = HoltenWaterConsts.a
    b = HoltenWaterConsts.b
    d = HoltenWaterConsts.d

    Tc = 228.2
    p0 = -300e6
    ρ0 = 1081.6482
    Rm = 461.523087
    τ = T/Tc
    _π = (p - p0)/(Rm*Tc*ρ0)
    lnπ,lnτ = log(_π), log(τ)
    ĝᴬ = zero(Base.promote_eltype(lnπ,lnτ))

    for k in eachindex(c)
        ĝᴬ += c[k]*exp(lnπ*b[k] + lnτ*a[k] - _π*d[k])
    end

    return ĝᴬ
end

p_scale(model::HoltenWaterModel,z) = 101325.0
T_scale(model::HoltenWaterModel,z) = 228.2

const HoltenWaterConsts = (
    c = [-8.1570681381655, 1.2875032, 7.0901673598012, -0.032779161, 0.73703949, -0.21628622, -5.1782479, 0.00042293517, 0.023592109, 4.3773754, -0.002996777, -0.96558018, 3.7595286, 1.2632441, 0.28542697, -0.85994947, -0.32916153, 0.090019616, 0.081149726, -3.2788213],
    a = [0.0, 0.0, 1.0, -0.2555, 1.5762, 1.64, 3.6385, -0.3828, 1.6219, 4.3287, 3.4763, 5.1556, -0.3593, 5.0361, 2.9786, 6.2373, 4.046, 5.3558, 9.0157, 1.2194],
    b = [0.0, 1.0, 0.0, 2.1051, 1.1422, 0.951, 0.0, 3.6402, 2.076, -0.0016, 2.2769, 0.0008, 0.3706, -0.3975, 2.973, -0.318, 2.9805, 2.9265, 0.4456, 0.1298],
    d = [0.0, 0.0, 0.0, -0.0016, 0.6894, 0.013, 0.0002, 0.0435, 0.05, 0.0004, 0.0528, 0.0147, 0.8584, 0.9924, 1.0041, 1.0961, 1.0228, 1.0303, 1.618, 0.5213],
)

export HoltenWater
