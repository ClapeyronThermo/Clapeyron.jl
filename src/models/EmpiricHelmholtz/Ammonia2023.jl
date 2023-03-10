struct Ammonia2023Consts <: EoSParam 
    Mw::Float64
    T_c::Float64
    P_c::Float64
    rho_c::Float64
    v_c::Float64
    Ttp::Float64
    ptp::Float64
    rhov_tp::Float64
    rhol_tp::Float64
    n::Vector{Float64}
    t::Vector{Float64}
    d::Vector{Int}
    l::Vector{Int}
    beta::Vector{Float64}
    gamma::Vector{Float64}
    eta::Vector{Float64}
    epsilon::Vector{Float64}
    bb::Vector{Float64}
    function Ammonia2023Consts()
        Mw = 17.03052 #g·mol-1
        T_c = 405.56    #K
        P_c = 11.3634e6 #Pa
        rho_c= 13696.0 # mol·m-3
        v_c = 1/rho_c
        Ttp = 195.49 #K
        ptp = 6.05339e3 
        rhov_tp  = 0.003740e3 
        rhol_tp = 43.091e3

        n = [0.006132232,1.7395866,-2.2261792,-0.30127553,0.08967023,-0.076387037,-0.84063963,-0.27026327,6.212578,-5.7844357,2.4817542,-2.3739168,0.01493697,-3.7749264,0.0006254348,
        -1.7359e-05,-0.13462033,0.07749072839,-1.6909858,0.93739074]
        t = [1.0,0.382,1.0,1.0,0.677,2.915,3.51,1.063,0.655,1.3,3.1,1.4395,1.623,0.643,1.13,4.5,1.0,4.0,4.3315,4.015]
        d = [4,1,1,2,3,3,2,3,1,1,1,2,2,1,3,3,1,1,1,1]
        l = [0,0,0,0,0,2,2,1]
        η = vcat(zeros(8),[0.42776,0.6424,0.8175,0.7995,0.91,0.3574,1.21,4.14,22.56,22.68,2.8452,2.8342])
        β = vcat(zeros(8),[1.708,1.4865,2.0915,2.43,0.488,1.1,0.85,1.14,945.64,993.85,0.3696,0.2962])
        γ = vcat(zeros(8),[1.036,1.2777,1.083,1.2906,0.928,0.934,0.919,1.852,1.05897,1.05277,1.108,1.313])
        ε = vcat(zeros(8),[-0.0726,-0.1274,0.7527,0.57,2.2,-0.243,2.96,3.02,0.9574,0.9576,0.4478,0.44689])
        b = vcat(zeros(18),[1.244,0.6826])

        return new(Mw,T_c,P_c,rho_c,v_c,
        Ttp,ptp,rhov_tp,rhol_tp,
        n,t,d,l,β,γ,η,ε,b)
    end
end

struct Ammonia2023 <: EmpiricHelmholtzModel
    components::Vector{String}
    consts::Ammonia2023Consts
    references::Vector{String}
end

"""
    Ammonia2023 <: EmpiricHelmholtzModel
    Ammonia2023()

## Input parameters

None

## Description
Ammonia Reference Equation of State (2023)
```
δ = ρ/ρc
τ = T/Tc
a⁰(δ,τ) = log(δ) + n⁰₁ + n⁰₂τ + n⁰₃log(τ) + ∑n⁰ᵢ(1-exp(-γ⁰ᵢτ)), i ∈ 4:7
aʳ(δ,τ)  = aʳ₁+ aʳ₂ + aʳ₃
aʳ₁(δ,τ)  =  ∑nᵢδ^(dᵢ)τ^(tᵢ), i ∈ 1:5
aʳ₂(δ,τ)  =  ∑nᵢexp(-δ^cᵢ)δ^(dᵢ)τ^(tᵢ), i ∈ 6:8
aʳ₃(δ,τ)  =  ∑nᵢexp(-ηᵢ(δ - εᵢ)^2 - βᵢ(τ - γᵢ)^2)δ^(dᵢ)τ^(tᵢ), i ∈ 9:18
aʳ₃(δ,τ)  =  ∑nᵢexp(-ηᵢ(δ - εᵢ)^2 - 1/(βᵢ*(τ -γᵢ)^2 + bᵢ))δ^(dᵢ)τ^(tᵢ), i ∈ 19:20

```
parameters  `n⁰`,`γ⁰`,`n`,`t`,`d`,`c`,`η`,`β`,`γ`,`ε` where obtained via fitting.

## References
1. Gao, K., Wu, J., Bell, I. H., Harvey, A. H., & Lemmon, E. W. (2023). A reference equation of state with an associating term for the thermodynamic properties of ammonia. Journal of Physical and Chemical Reference Data, 52(1), 013102. [doi:10.1063/5.0128269](https://doi.org/10.1063/5.0128269)
"""
Ammonia2023

Ammonia2023() = Ammonia2023(["ammonia"],Ammonia2023Consts(),["10.1063/5.0128269"])


is_splittable(::Ammonia2023) = false

function _f0(::Ammonia2023,δ,τ)
    _1 = one(τ)
    a₁ = -6.59406093943886
    a₂ = 5.601011519879   
    u = (1646.0,3965.0,7231.0)
    v = (2.224,3.148,0.9579)
    α₀ = log(δ) + 3*log(τ) + a₁ + a₂*τ +
     v[1]*log(_1-exp(-u[1]*τ)) + 
     v[2]*log(_1-exp(-u[2]*τ)) + 
     v[3]*log(_1-exp(-u[3]*τ))
     return α₀
end

function _fr1(model::Ammonia2023,δ,τ)
    n = model.consts.n
    t = model.consts.t
    d = model.consts.d
    l = model.consts.l
    η = model.consts.eta
    β = model.consts.beta
    γ = model.consts.gamma
    ε = model.consts.epsilon
    b = model.consts.bb

    αᵣ = zero(δ+τ)
    
    lnδ = log(δ)
    lnτ = log(τ)
    for k in 1:5
        αᵣ += n[k]*exp(lnδ*d[k] + lnτ*t[k])
    end

    for k in 6:8
        αᵣ += n[k]*exp(lnδ*d[k] + lnτ*t[k] - δ^l[k])
    end

    for k in 9:18
        αᵣ += n[k]*exp(lnδ*d[k] + lnτ*t[k] - η[k]*(δ-ε[k])^2  - β[k]*(τ-γ[k])^2)
    end
  
    for k in 19:20
        αᵣ += n[k]*exp(lnδ*d[k] + lnτ*t[k] - η[k]*(δ-ε[k])^2  + 1/(β[k]*(τ-γ[k])^2 + b[k]))
    end

    return αᵣ
end
#ancillary equations for calculation of P_sat, T_sat rhovsat y rholsat
function _Ammonia2023_psat(T)
    T_c = 405.56
    P_c = 11.3634e6
    T>T_c && return zero(T)/zero(T)
    Tr = T/T_c
    θ = 1.0-Tr
    lnPsatPc = evalexppoly(θ,(-7.3128,3.8888,-2.9908,-2.8636),(1.0,1.5,1.6,1.7))
    Psat = exp(lnPsatPc)*P_c
    return Psat
end

function _Ammonia2023_tsat(p)
    P_c = 11.3634e6
    T_c = 405.56
    p > P_c && return zero(p)/zero(p)
    #first aproximation

    if p < 1e5 #stull, 1947 antoine coeff for ammonia
        A,B,C = 3.18757,506.713,-80.78
    else
        A,B,C = 4.86886,1113.928,-10.409
    end

    T0 = B/(A - log10(p*1e-5)) - C
    T0 > T_c && (T0 = T_c*p/P_c)
    f(T) = _Ammonia2023_psat(T) - p
    prob = Roots.ZeroProblem(f,T0)
    return Roots.solve(prob,Roots.Order0())
end

function _Ammonia2023_rholsat(T)
    T_c = 405.56
    ρ_c =13696.0
    T>T_c && return zero(T)/zero(T)
    Tr = T/T_c
    θ = 1.0-Tr
    ρ_r = evalexppoly(θ,(0.051236,3.7925,-3.5929,4.6409,-1.9893,1.5978),(0.07,0.46,0.77,1.05,1.25,8.0)) + 1.0
    ρ_l = ρ_r*ρ_c
    return ρ_l
end

#(-0.089966,-3.8722,-8.1183,-25.293,-54.279,-400.83)

function _Ammonia2023_rhovsat(T)
    T_c = 405.56
    ρ_c =13696.0
    T>T_c && return zero(T)/zero(T)
    Tr = T/T_c
    θ = 1.0 - Tr
    log_ρ_v_ρ_c =evalexppoly(θ,(-0.089966,-3.8722,-8.1183,-25.293,-54.279,-400.83),(0.112,0.473,1.5,3.875,8.0,20.0))
    ρ_v = exp(log_ρ_v_ρ_c)*ρ_c
    return ρ_v
end

function a_ideal(model::Ammonia2023,V,T,z=SA[1.])
    T_c = model.consts.T_c
    rho_c = model.consts.rho_c
    N = only(z)
    rho = (N/V)
    δ = rho/rho_c
    τ = T_c/T
    return  _f0(model,δ,τ)
end

function a_res(model::Ammonia2023,V,T,z=SA[1.])
    T_c = model.consts.T_c
    rho_c = model.consts.rho_c
    N = only(z)
    rho = (N/V)
    δ = rho/rho_c
    τ = T_c/T
    return  _fr1(model,δ,τ)
end

function eos(model::Ammonia2023, V, T, z=SA[1.0];phase=:unknown)
    R = 8.314462618
    T_c = model.consts.T_c
    rho_c = model.consts.rho_c
    N = only(z)
    rho = (N/V)
    δ = rho/rho_c
    τ = T_c/T
    return N*R*T*(_f0(model,δ,τ)+_fr1(model,δ,τ))
end

function eos_res(model::Ammonia2023,V,T,z=SA[1.0];phase=:unknown)
    R = 8.314462618
    T_c = model.consts.T_c
    rho_c = model.consts.rho_c
    N = only(z)
    rho = (N/V)
    δ = rho/rho_c
    τ = T_c/T
    return N*R*T*_fr1(model,δ,τ)
end

mw(model::Ammonia2023) = SA[model.consts.Mw]

molecular_weight(model::Ammonia2023,z = @SVector [1.]) = model.consts.Mw*0.001

T_scale(model::Ammonia2023,z=SA[1.0]) = model.consts.T_c

p_scale(model::Ammonia2023,z=SA[1.0]) = model.consts.P_c

lb_volume(model::Ammonia2023,z=SA[1.0]) = 1/53130 #finally, an eos model that mentions it max density.

Base.length(::Ammonia2023) = 1

function Base.show(io::IO,mime::MIME"text/plain",model::Ammonia2023)
    print(io,"Ammonia Reference Equation of State")
end

function x0_sat_pure(model::Ammonia2023,T,z=SA[1.0])
    vv = 1.0/_Ammonia2023_rhovsat(T)
    vl = 1.0/_Ammonia2023_rholsat(T)
    return (vl,vv)
end

function x0_volume_liquid(model::Ammonia2023,T,z = SA[1.0])
    return  1/_Ammonia2023_rholsat(min(T,model.consts.T_c*one(T)))
end

x0_psat(model::Ammonia2023,T,crit=nothing) = _Ammonia2023_psat(T)

function x0_saturation_temperature(model::Ammonia2023,p)
    T = _Ammonia2023_tsat(p)
    vl = 1.0/_Ammonia2023_rholsat(T)
    vv = 1.0/_Ammonia2023_rhovsat(T)
    return (T,vl,vv)
end

export Ammonia2023