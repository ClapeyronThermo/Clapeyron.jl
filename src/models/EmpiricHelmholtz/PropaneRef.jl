struct PropaneRefConsts <: EoSParam 
    Mw::Float64
    T_c::Float64
    P_c::Float64
    rho_c::Float64
    v_c::Float64
    n::Vector{Float64}
    t::Vector{Float64}
    d::Vector{Int}
    l::Vector{Int}
    beta::Vector{Float64}
    gamma::Vector{Float64}
    eta::Vector{Float64}
    epsilon::Vector{Float64}
    function PropaneRefConsts()
        Mw = 44.09562 #g·mol-1
        T_c = 369.89    #K
        P_c = 4.2512e6 #Pa
        rho_c= 5000.0 # mol·m-3
        v_c = 1/rho_c
        n = [0.042910051, 1.7313671, -2.4516524, 0.34157466, -0.46047898, -0.66847295, 0.20889705, 0.19421381, -0.22917851, -0.60405866, 0.066680654, 0.017534618, 0.33874242, 0.22228777, -0.23219062, -0.09220694, -0.47575718, -0.017486824]
        t = [1.0, 0.33, 0.8, 0.43, 0.9, 2.46, 2.09, 0.88, 1.09, 3.25, 4.62, 0.76, 2.5, 2.75, 3.05, 2.55, 8.4, 6.75]
        d = [4,1,1,2,2,1,3,6,6,2,3,1,1,1,2,2,4,1]
        l = [1,1,1,1,2,2]
        η = [0.963,1.977,1.917,2.307,2.546,3.28,14.6]
        β = [2.33,3.47,3.15,3.19,0.92,18.8,547.8]
        γ = [0.684,0.829,1.419,0.817,1.500,1.426,1.093]
        ε = [1.283,0.6936,0.788,0.473,0.8577,0.271,0.948]
        return new(Mw,T_c,P_c,rho_c,v_c,
        n,t,d,l,β,γ,η,ε)
    end
end

struct PropaneRef <: EmpiricHelmholtzModel
    components::Vector{String}
    consts::PropaneRefConsts
    references::Vector{String}
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
aʳ₁(δ,τ)  =  ∑nᵢδ^(dᵢ)τ^(tᵢ), i ∈ 1:5
aʳ₂(δ,τ)  =  ∑nᵢexp(-δ^cᵢ)δ^(dᵢ)τ^(tᵢ), i ∈ 6:11
aʳ₃(δ,τ)  =  ∑nᵢexp(-ηᵢ(δ - εᵢ)^2 - βᵢ(τ - γᵢ)^2)δ^(dᵢ)τ^(tᵢ), i ∈ 12:18

```
parameters  `n⁰`,`γ⁰`,`n`,`t`,`d`,`c`,`η`,`β`,`γ`,`ε` where obtained via fitting.

## References
1. Lemmon, E. W., McLinden, M. O., & Wagner, W. (2009). Thermodynamic properties of propane. III. A reference equation of state for temperatures from the melting line to 650 K and pressures up to 1000 MPa. Journal of Chemical and Engineering Data, 54(12), 3141–3180. [doi:10.1021/je900217v](https://doi.org/10.1021/je900217v)
"""
PropaneRef

PropaneRef() = PropaneRef(["propane"],PropaneRefConsts(),["1021/je900217v"])


is_splittable(::PropaneRef) = false

function _f0(::PropaneRef,δ,τ)
    δ,τ = promote(δ,τ)
    _1 = one(δ)
    n₁ = -4.970583
    n₂ = 4.29352     
    γ = (1.062478, 3.344237,5.363757,11.762957)
    n = (3.043,5.874,9.337,7.922)
    α₀ = log(δ) + 3*log(τ) + n₁ + n₂*τ +
     n[1]*log(_1-exp(-γ[1]*τ)) + 
     n[2]*log(_1-exp(-γ[2]*τ)) + 
     n[3]*log(_1-exp(-γ[3]*τ)) + 
     n[4]*log(_1-exp(-γ[4]*τ))
     return α₀
end


function _fr1(model::PropaneRef,δ,τ)
    δ,τ = promote(δ,τ)
    n = model.consts.n
    t = model.consts.t
    d = model.consts.d
    l = model.consts.l
    η = model.consts.eta
    β = model.consts.beta
    γ = model.consts.gamma
    ε = model.consts.epsilon
    αᵣ = zero(δ)
    
    for k in 1:5
        αᵣ += n[k]*(δ^d[k])*(τ^t[k])
    end

    for k in 6:11
        αᵣ += n[k]*(δ^d[k])*(τ^t[k])*exp(-δ^l[k-5])
    end

    for k in 12:18
        _k = k - 11
        mul = exp(-η[_k]*(δ-ε[_k])^2  - β[_k]*(τ-γ[_k])^2)
        αᵣ += n[k]*(δ^d[k])*(τ^t[k])*mul
    end
    return αᵣ
end
#ancillary equations for calculation of P_sat, T_sat rhovsat y rholsat
function _propaneref_psat(T)
    T_c = 369.89
    P_c = 4.2512e6
    T>T_c && return zero(T)/zero(T)
    Tr = T/T_c
    θ = 1.0-Tr
    lnPsatPc = (-6.7722*θ + 1.6938*θ^1.5 -1.3341*θ^2.2 -3.1876*θ^4.8 + 0.94937*θ^6.2)/Tr
    Psat = exp(lnPsatPc)*P_c
    return Psat
end

function _propaneref_tsat(p)
    P_c = 4.2512e6
    p > P_c && return zero(p)/zero(p)
    #first aproximation
    A,B,C = 13.6515,1850.8,249.99-273.15
    T0 = B/(A - log(p*1e-3)) - C
    T0 > 369.89 && (T0 = 369.89*p/P_c)
    f(T) = _propaneref_psat(T) - p
    prob = Roots.ZeroProblem(f,T0)
    return Roots.solve(prob,Roots.Order0())
end

function _propaneref_rholsat(T)
    T_c = 369.89
    ρ_c = 5000.0
    T>T_c && return zero(T)/zero(T)
    Tr = T/T_c
    θ = 1.0-Tr
    ρ_l = (1.0 + 1.82205*θ^0.345 + 0.65802*θ^0.74 + 0.21109*θ^2.6 + 0.083973*θ^7.2)*ρ_c
    return ρ_l
end

function _propaneref_rhovsat(T)
    T_c = 369.89
    ρ_c = 5000.0
    T>T_c && return zero(T)/zero(T)
    Tr = T/T_c
    θ = 1.0 - Tr
    log_ρ_v_ρ_c = (-2.4887*θ^0.3785 -5.1069*θ^1.07 -12.174*θ^2.7 -30.495*θ^5.5 -52.192*θ^10 -134.89*θ^20)
    ρ_v = exp(log_ρ_v_ρ_c)*ρ_c
    return ρ_v
end

function a_ideal(model::PropaneRef,V,T,z=SA[1.])
    T_c = model.consts.T_c
    rho_c = model.consts.rho_c
    N = only(z)
    rho = (N/V)
    δ = rho/rho_c
    τ = T_c/T
    return  _f0(model,δ,τ)
end

function a_res(model::PropaneRef,V,T,z=SA[1.])
    T_c = model.consts.T_c
    rho_c = model.consts.rho_c
    N = only(z)
    rho = (N/V)
    δ = rho/rho_c
    τ = T_c/T
    return  _fr1(model,δ,τ)
end

function eos(model::PropaneRef, V, T, z=SA[1.0];phase=:unknown)
    R =8.314472
    T_c = model.consts.T_c
    rho_c = model.consts.rho_c
    N = only(z)
    rho = (N/V)
    δ = rho/rho_c
    τ = T_c/T
    return N*R*T*(_f0(model,δ,τ)+_fr1(model,δ,τ))
end

function eos_res(model::PropaneRef,V,T,z=SA[1.0];phase=:unknown)
    R =8.314472
    T_c = model.consts.T_c
    rho_c = model.consts.rho_c
    N = only(z)
    rho = (N/V)
    δ = rho/rho_c
    τ = T_c/T
    return N*R*T*_fr1(model,δ,τ)
end

mw(model::PropaneRef) = SA[model.consts.Mw]

molecular_weight(model::PropaneRef,z = @SVector [1.]) = model.consts.Mw*0.001

T_scale(model::PropaneRef,z=SA[1.0]) = model.consts.T_c

p_scale(model::PropaneRef,z=SA[1.0]) = model.consts.P_c

lb_volume(model::PropaneRef,z=SA[1.0]) = 6.0647250138479785e-5 #calculated at 1000 MPa and 650 K

Base.length(::PropaneRef) = 1

function Base.show(io::IO,mime::MIME"text/plain",model::PropaneRef)
    print(io,"Propane Reference Equation of State")
end

function x0_sat_pure(model::PropaneRef,T,z=SA[1.0])
    vv = 1.0/_propaneref_rhovsat(T)
    vl = 1.0/_propaneref_rholsat(T)
    return (vl,vv)
end

function x0_volume_liquid(model::PropaneRef,T,z = SA[1.0])
    return  1/_propaneref_rholsat(min(T,369.88889*one(T)))
end

x0_psat(model::PropaneRef,T,crit=nothing) = _propaneref_psat(T)

function x0_saturation_temperature(model::PropaneRef,p)
    T = _propaneref_tsat(p)
    vl = 1.0/_propaneref_rholsat(T)
    vv = 1.0/_propaneref_rhovsat(T)
    return (T,vl,vv)
end

export PropaneRef

#=
    R = 8.314472 #J·mol-1·K-1
    ,Mw = 44.09562 #g·mol-1
    ,T_c = 369.89    #K
    ,P_c = 4.2512e6 #Pa
    ,rho_c= 5000.0 # mol·m-3
    ,T_triple = 85.525 #K
    ,P_Triple = 0.00017 #Pa
    ,rho_v_triple = 2.4e-07 # mol·m-3
    ,rho_l_triple = 16626.0 # mol·m-3
    ,T_0 = 231.036 #K
    ,rho_v_0 = 54.8 # mol·m-3
    ,rho_l_0 = 13173.0 # mol·m-3
    ,ref_T_ideal = 273.15
    ,ref_P_0 = 100000.0
    ,ref_T_ideal_H = 26148.48 ## J·mol-1
    ,ref_T_ideal_S = 157.9105 ##J·mol-1·K-1
=#