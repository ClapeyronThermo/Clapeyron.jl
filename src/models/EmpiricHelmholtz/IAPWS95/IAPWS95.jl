struct IAPWS95Consts <: EoSParam
    Tc::Float64
    Pc::Float64
    Vc::Float64
    Mw::Float64
    n::Vector{Float64}
    t::Vector{Float64}
    d::Vector{Int}
    c::Vector{Int}
    function IAPWS95Consts()
        n = [0.12533547935523e-1, 0.78957634722828e1, -0.87803203303561e1,
        0.31802509345418, -0.26145533859358, -0.78199751687981e-2,
        0.88089493102134e-2,

        -0.66856572307965, 0.20433810950965, -0.66212605039687e-4,
        -0.19232721156002, -0.25709043003438, 0.16074868486251,
        -0.040092828925807, 0.39343422603254e-6, -0.75941377088144e-5,
        0.56250979351888e-3, -0.15608652257135e-4, 0.11537996422951e-8,
        .36582165144204e-6, -.13251180074668e-11, -.62639586912454e-9,
        -0.10793600908932, 0.17611491008752e-1, 0.22132295167546,
        -0.40247669763528, 0.58083399985759, 0.49969146990806e-2,
        -0.31358700712549e-1, -0.74315929710341, 0.47807329915480,
        0.20527940895948e-1, -0.13636435110343, 0.14180634400617e-1,
        0.83326504880713e-2, -0.29052336009585e-1, 0.38615085574206e-1,
        -0.20393486513704e-1, -0.16554050063734e-2, .19955571979541e-2,
        0.15870308324157e-3, -0.16388568342530e-4, 0.43613615723811e-1,
        0.34994005463765e-1, -0.76788197844621e-1, 0.22446277332006e-1,
        -0.62689710414685e-4, -0.55711118565645e-9, -0.19905718354408,
        0.31777497330738, -0.11841182425981]

        d = [1, 1, 1, 2, 2, 3, 4,

        1, 1, 1, 2, 2, 3, 4, 4, 5, 7, 9, 10, 11, 13, 15, 1, 2, 2, 2, 3,
        4, 4, 4, 5, 6, 6, 7, 9, 9, 9, 9, 9, 10, 10, 12, 3, 4, 4, 5, 14, 3, 6, 6, 6]

        t = [-0.5, 0.875, 1, 0.5, 0.75, 0.375, 1,
        4, 6, 12, 1, 5, 4, 2, 13, 9, 3, 4, 11, 4, 13, 1, 7, 1, 9, 10, 10, 3, 7,
        10, 10, 6, 10, 10, 1, 2, 3, 4, 8, 6, 9, 8, 16, 22, 23,23, 10, 50, 44, 46, 50]

        c = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 6, 6,6, 6]
        return new(647.096,2.2064e7,5.594803726708074e-5,18.015268,
        n,t,d,c)
    end
end

struct IAPWS95 <: EmpiricHelmholtzModel
    components::Vector{String}
    consts::IAPWS95Consts
    references::Vector{String}
end

function crit_pure(model::IAPWS95)
    return (model.consts.Tc,model.consts.Pc,model.consts.Vc)
end



"""
    IAPWS95 <: EmpiricHelmholtzModel
    IAPWS95()

## Input parameters

None

## Description

IAPWS95 (International Association for the Properties of Water and Steam) Pure water Model, 2018 update.

```
δ = ρ/ρc
τ = T/Tc
a⁰(δ,τ) = log(δ) + n⁰₁ + n⁰₂τ + n⁰₃log(τ) + ∑n⁰ᵢ(1-exp(-γ⁰ᵢτ)), i ∈ 4:8
aʳ(δ,τ)  = aʳ₁+ aʳ₂ + aʳ₃ + aʳ₄
aʳ₁(δ,τ)  =  ∑nᵢδ^(dᵢ)τ^(tᵢ), i ∈ 1:7
aʳ₂(δ,τ)  =  ∑nᵢexp(-δ^cᵢ)δ^(dᵢ)τ^(tᵢ), i ∈ 8:51
aʳ₃(δ,τ)  =  ∑nᵢexp(-αᵢ(δ - εᵢ)^2 - βᵢ(τ - γᵢ)^2)δ^(dᵢ)τ^(tᵢ), i ∈ 52:54
aʳ₄(δ,τ) = ∑nᵢδΨΔ^(bᵢ), i ∈ 55:56
Δ = θ^2 + Bᵢ[(δ - 1)^2]^aᵢ
θ = (1 - τ) + Aᵢ[(δ - 1)^2]^(1/2βᵢ)
Ψ = exp(-Cᵢ(δ - 1)^2 - Dᵢ(τ - 1)^2)
```
parameters `n⁰`,`γ⁰`,`n`,`t`,`d`,`c`,`α`,`β`,`γ`,`ε`,`A`,`B`,`C`,`D` where obtained via fitting.

## References

1. Wagner, W., & Pruß, A. (2002). The IAPWS formulation 1995 for the thermodynamic properties of ordinary water substance for general and scientific use. Journal of physical and chemical reference data, 31(2), 387–535. [doi:10.1063/1.1461829](https://doi.org/10.1063/1.1461829)
2. IAPWS R6-95 (2018). Revised Release on the IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary Water Substance for General and Scientific Use

"""
IAPWS95

IAPWS95() = IAPWS95(["water"],IAPWS95Consts(),["IAPWS R6-95(2018)"])

function iapws_f0(model,δ,τ)
    n= (-8.3204464837497, 6.6832105275932, 3.00632,0.012436, 0.97315, 1.2795, 0.96956, 0.24873)
    γ = (0.0, 0.0, 0.0, 1.28728967, 3.53734222, 7.74073708, 9.24437796,27.5075105)
    res = log(δ)+n[1] + n[2]*τ + n[3]*log(τ)
    @inbounds for i = 4:8
        res += n[i]*log(-expm1(-γ[i]*τ))
    end
    return res
end

_f0(model::IAPWS95,δ,τ) = iapws_f0(model,δ,τ)

function _fr(model::IAPWS95,δ,τ)
    n = model.consts.n::Vector{Float64}
    d = model.consts.d::Vector{Int64}
    t = model.consts.t::Vector{Float64}
    c = model.consts.c::Vector{Int64}
    res=zero(δ+τ)
    for i = 1:7
        res += n[i]* (δ^d[i]) * (τ^t[i])
    end

    for i = 8:51
        ic = i-7
        res += n[i]* (δ^d[i]) * (τ^t[i]) * exp(-δ^c[ic])
    end
    
    nτt1,nτt2,nτt3 = (-0.31306260323435e2, 0.31546140237781e2*τ, -0.25213154341695e4*τ^4)
    δd = δ^3
    #for i = 52:54
    #(δ-ε)^2
    _δ = (δ-1.0)^2#
    gauss_d = 20*_δ
    gauss_1 = exp(-150*((τ-1.21)^2)  -gauss_d)
    gauss_2 = gauss_1
    gauss_3 = exp(-250*((τ-1.25)^2)  -gauss_d)

    res += δd*(nτt1*gauss_1 + nτt2*gauss_2 + nτt3*gauss_3)

    #for i = 55:56

    nδ1,nδ2 = (-0.14874640856724*δ, 0.31806110878444*δ)
    _τ = (τ-1.0)^2
    Θ = (1.0-τ) + 0.32*_δ^(1.6666666666666667)
    Δ = Θ^2 + 0.2*_δ^3.5
    Ψ1 = exp(- 28*_δ - 700*_τ)
    Ψ2 = exp(- 32*_δ - 800*_τ)
    Δb1 = Δ^0.85
    Δb2 = Δ^0.95
    res += nδ1*Δb1*Ψ1 + nδ2*Δb2*Ψ2

    return res

end

function a_ideal(model::IAPWS95,V,T,z=SA[1.0])
    Σz = only(z) #single component
    v = V/Σz
     mass_v =  v*1000.0*0.055508472036052976
     rho = one(mass_v)/mass_v
     δ = rho*0.003105590062111801 #/322
     τ = 647.096/T
     return 0.9999890238768239*_f0(model,δ,τ)
end

function a_res(model::IAPWS95,V,T,z=SA[1.0])
    Σz = only(z) #single component
    v = V/Σz
    #R value calculated from molecular weight and specific gas constant
     #return 8.3143713575874*T*_f(model, molar_to_weight(1/v,[model.molecularWeight],[1.0]),T)
     #println(molar_to_weight(1/v,[model.molecularWeight],[1.0]))'
     mass_v =  v*1000.0*0.055508472036052976
     rho = one(mass_v)/mass_v
     δ = rho*0.003105590062111801 #/322
     τ = 647.096/T
     return 0.9999890238768239*_fr(model,δ,τ)
end


function eos(model::IAPWS95, V, T, z=SA[1.0])
    Σz = only(z) #single component
    v = V/Σz
    #R value calculated from molecular weight and specific gas constant
     #return 8.3143713575874*T*_f(model, molar_to_weight(1/v,[model.molecularWeight],[1.0]),T)
     #println(molar_to_weight(1/v,[model.molecularWeight],[1.0]))'
     mass_v =  v*1000.0*0.055508472036052976
     rho = one(mass_v)/mass_v
     δ = rho*0.003105590062111801 #/322
     τ = 647.096/T
    return R̄*Σz*T*0.9999890238768239*(_fr(model,δ,τ)+_f0(model,δ,τ))
end

function eos_res(model::IAPWS95, V, T, z=SA[1.0])
    Σz = only(z) #single component
    v = V/Σz
    #R value calculated from molecular weight and specific gas constant
     #return 8.3143713575874*T*_f(model, molar_to_weight(1/v,[model.molecularWeight],[1.0]),T)
     #println(molar_to_weight(1/v,[model.molecularWeight],[1.0]))'
     mass_v =  v*1000.0*0.055508472036052976
     rho = one(mass_v)/mass_v
     δ = rho*0.003105590062111801 #/322
     τ = 647.096/T
    return R̄*Σz*T*0.9999890238768239*_fr(model,δ,τ)
end

const watersat_data = (;n = [0.116_705_214_527_67E4,
    -0.724_213_167_032_06E6,
    -0.170_738_469_400_92E2,
    0.120_208_247_024_70E5,
    -0.323_255_503_223_33E7,
    0.149_151_086_135_30E2,
    -0.482_326_573_615_91E4,
    0.405_113_405_420_57E6,
    -0.238_555_575_678_49,
    0.650_175_348_447_98E3]
)

#psat and tsat are modified from SteamTables.jl
function water_p_sat(t)
    t > 647.096 && return zero(t)/zero(t)
    n = watersat_data.n
    Θ = t + n[9]/(t - n[10])
    A =      Θ^2 + n[1]*Θ + n[2]
    B = n[3]*Θ^2 + n[4]*Θ + n[5]
    C = n[6]*Θ^2 + n[7]*Θ + n[8]
    P = (2C / (-B + √(B^2 - 4*A*C)))^4
    return P*1000000
end

function water_t_sat(p)
    p > 2.2064e7 && return zero(p)/zero(p)
    n = watersat_data.n
    P = p/1000000
    β = P^0.25
    E =      β^2 + n[3]*β + n[6]
    F = n[1]*β^2 + n[4]*β + n[7]
    G = n[2]*β^2 + n[5]*β + n[8]
    D = 2G / (-F - √(F^2 - 4*E*G))
    T = (n[10]+D-√((n[10]+D)^2 - 4(n[9]+n[10]*D)))/2.0
    return T
end

function x0_saturation_temperature(model::IAPWS95,p)
    T = water_t_sat(p)
    vl = saturated_water_liquid(T)
    vv = saturated_water_vapor(T)
    return (T,vl,vv)
end

x0_psat(model::IAPWS95,T,crit=nothing) = water_p_sat(T)

#from MoistAir.jl, liquid
function saturated_water_liquid(Tk)
    ρm= ( -0.2403360201e4 +
        Tk*(-0.140758895e1 +
        Tk * (0.1068287657e0 +
        Tk*(-0.2914492351e-3 + Tk*(0.373497936e-6 - Tk*0.21203787e-9))))) /
        ( -0.3424442728e1 + 0.1619785e-1*Tk )
    #kg/m3 -> m3/kg *(kg/mol)
    mw = 18.015268*0.001
    v = mw/ρm
end

function saturated_water_vapor(Tk)
    Blin = 0.70e-8 - 0.147184e-8 * exp(1734.29/Tk)
    Clin = 0.104e-14 - 0.335297e-17*exp(3645.09/Tk)
    p = water_p_sat(Tk)
    return R̄*Tk/p * (1.0 + p*(Blin  + p*Clin))
end

mw(model::IAPWS95) = model.consts.Mw
molecular_weight(model::IAPWS95,z=SA[1.]) = model.consts.Mw*1e-3

function x0_volume(model::IAPWS95,p,T,z=[1.0];phase = :unknown)
    if phase == :unknown || is_liquid(phase)
        if model.consts.Pc > p
            return sat_v = saturated_water_liquid(T)
        else
            return model.consts.Vc #must look for better initial point here
        end
    elseif is_vapour(phase)
        return 1.1*saturated_water_vapor(T)
    elseif is_supercritical(phase)
        return model.consts.Vc
    else
        _0 = zero(p+T+first(z))
        return _0/_0
    end
end

function x0_sat_pure(model::IAPWS95,T)
    if T > 647.096
        nan = zero(T)/zero(T)
        return nan,nan
    else
        vl = saturated_water_liquid(T)
        vg = saturated_water_vapor(T)
        return (vl,vg)
    end
end

function T_scale(model::IAPWS95,z=SA[1.0])
    return model.consts.Tc
end

function p_scale(model::IAPWS95,z=SA[1.0])
    return model.consts.Pc
end

Base.length(::IAPWS95) = 1

function Base.show(io::IO,model::IAPWS95)
    return eosshow(io,model)
end

function Base.show(io::IO,mime::MIME"text/plain",model::IAPWS95)
    return eosshow(io,mime,model)
end

lb_volume(model::IAPWS95, z=SA[1.0]; phase = :unknown) = 1.4393788065379039e-5

export IAPWS95,IAPWS95Ideal


