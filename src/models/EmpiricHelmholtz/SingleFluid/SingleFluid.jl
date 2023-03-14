struct EmpiricSingleFluid{𝔾,ℙ,𝔸,𝕀,ℝ} <: EmpiricHelmholtzModel
    type::𝔾
    components::Vector{String}
    properties::ℙ
    ancilliaries::𝔸
    ideal::𝕀
    residual::ℝ
    references::Vector{String}
end

## Description
"""
Single Multiparameter Fluid Equation of state.

```
δ = ρ/ρc
τ = T/Tc
a⁰(δ,τ)   =  log(δ) + a₁ + a₂τ + (c₀ - 1)*log(τ) + ∑vᵢ(1-exp(uᵢτ))
aʳ(δ,τ)   =  aʳ₁+ aʳ₂ + aʳ₃
aʳ₁(δ,τ)  =  ∑nᵢδ^(dᵢ)τ^(tᵢ), i ∈ k_pol
aʳ₂(δ,τ)  =  ∑nᵢexp(-δ^cᵢ)δ^(dᵢ)τ^(tᵢ), i ∈ k_exp
aʳ₃(δ,τ)  =  ∑nᵢexp(-ηᵢ(δ - εᵢ)^2 - βᵢ(τ - γᵢ)^2)δ^(dᵢ)τ^(tᵢ), i ∈ k_gauss
aʳ₃(δ,τ)  =  ∑nᵢexp(-ηᵢ(δ - εᵢ)^2 - 1/(βᵢ*(τ -γᵢ)^2 + bᵢ))δ^(dᵢ)τ^(tᵢ), i ∈ k_assoc
```

All parameters are fitted, to allow a equation of state of a single fluid with property calculations as close as possible to the experimental values.
"""

struct IdealEmpiricSingleFluid{𝔾,ℙ,𝕀} <: IdealModel
    type::𝔾
    components::Vector{String}
    properties::ℙ
    ideal::𝕀
    references::Vector{String}
end

function IdealEmpiricSingleFluid(model::EmpiricSingleFluid)
    return IdealEmpiricSingleFluid(model.type,model.components,model.properties,model.ideal,model.references)
end

struct EmpiricSingleFluidIdealParam <:EoSParam
    a1::Float64
    a2::Float64
    c0::Float64
    u::Vector{Float64}
    v::Vector{Float64}
end

struct EmpiricSingleFluidParam <: EoSParam
    iterators::Vector{UnitRange{Int}}
    n::Vector{Float64}
    t::Vector{Float64}
    d::Vector{Int}
    l::Vector{Int}
    eta::Vector{Float64}
    beta::Vector{Float64}
    gamma::Vector{Float64}
    epsilon::Vector{Float64}
    b_assoc::Vector{Float64}

    function EmpiricSingleFluidParam(n,t,d,l = Int[],eta = Float64[],beta = Float64[],gamma = Float64[], epsilon = Float64[],b_assoc = Float64[])
        param = new(Vector{UnitRange{Int}}(undef,0),n,t,d,l,eta,beta,gamma,epsilon,b_assoc)
        _calc_iterators!(param)
        return param
    end
end

function _calc_iterators!(param::EmpiricSingleFluidParam)
    n,t,d,l = param.n,param.t,param.d,param.l
    eta,beta,gamma,epsilon = param.eta,param.beta,param.gamma,param.epsilon
    b_assoc = param.b_assoc

    @assert length(n) == length(t) == length(d)
    @assert length(l) < length(d)
    @assert length(eta) == length(beta) == length(gamma) == length(epsilon)
    @assert length(b_assoc) < length(beta)

    #we start from the assoc term, backwards
    length_n = length(n)
    length_beta = length(beta)
    length_b = length(b_assoc)

    length_pol = length_n - length_beta - length(l)
    length_exp = length_n - length_beta
    length_gauss = length_n - length_b
    k_pol = 1:length_pol
    k_exp = (length_pol+1):length_exp
    k_gauss = (length_exp+1):length_gauss
    k_assoc = (length_gauss+1):length_n
    resize!(param.iterators,4)
    param.iterators .= (k_pol,k_exp,k_gauss,k_assoc)
    return param
end

function recombine_impl!(model::EmpiricSingleFluid)
    _calc_iterators!(model.residual)
    return model
end

struct EmpiricSingleFluidProperties <: EoSParam
    Mw::Float64 #Molecular Weight, g/mol
    Tc::Float64 #Critical temperature, K
    Pc::Float64 #Critical Pressure,Pa
    rhoc::Float64 #Critical density, mol/m3
    lb_volume::Float64 #lower bound volume, mol/m3
    Ttp::Float64 #triple point temperature, K
    ptp::Float64 #triple point pressure, Pa
    rhov_tp::Float64 #triple point vapor volume, mol/m3
    rhol_tp::Float64 #triple point liquid volume, mol/m3
    acentricfactor::Float64 #acentric factor
    Rgas::Float64 #gas constant used

    function EmpiricSingleFluidProperties(Mw,Tc,Pc,rhoc,lb_volume,
        Ttp = NaN,ptp = NaN, rhov_tp = NaN,rhol_tp = NaN, acentric_factor = NaN, Rgas = R̄)
        return new(Mw,Tc,Pc,rhoc,lb_volume, Ttp,ptp,rhov_tp,rhol_tp,acentric_factor,Rgas)
    end
end

idealmodel(model::EmpiricSingleFluid) = IdealEmpiricSingleFluid(model)

R_gas(model::EmpiricSingleFluid) = model.properties.Rgas
R_gas(model::IdealEmpiricSingleFluid) = model.properties.Rgas

function _f0(model::Union{EmpiricSingleFluid,IdealEmpiricSingleFluid},δ,τ)
    a₁ = model.ideal.a1
    a₂ = model.ideal.a2
    c₀ = model.ideal.c0

    α₀ = log(δ) + a₁ + a₂*τ + (c₀ - 1)*log(τ)
    

    u = model.ideal.u
    length(u) == 0 && return α₀
    v = model.ideal.v
    for i in eachindex(u)
        α₀ += v[i]*log(1 - exp(-u[i]*τ))
    end
    return α₀
end

function _fr1(model::EmpiricSingleFluid,δ,τ,type = model.type)

    αᵣ = zero(δ+τ)
    lnδ = log(δ)
    lnτ = log(τ)

    ℙ = model.residual

    n,t,d,l,η,β,γ,ε,b = ℙ.n,ℙ.t,ℙ.d,ℙ.l,ℙ.eta,ℙ.beta,ℙ.gamma,ℙ.epsilon,ℙ.b_assoc

    k_pol,k_exp,k_gauss,k_assoc = model.residual.iterators

    #strategy for storing.
    #n, t, d, gauss values, always require views
    #l, b does not require views. they are used just once.

    #Polynomial terms
    n_pol = view(n,k_pol)
    t_pol = view(t,k_pol)
    d_pol = view(d,k_pol)
    αᵣ += _fr1_pol(δ,τ,lnδ,lnτ,αᵣ,n_pol,t_pol,d_pol)

    #Exponential terms
    length(k_exp) == 0 && return αᵣ
    n_exp = view(n,k_exp)
    t_exp = view(t,k_exp)
    d_exp = view(d,k_exp)
    αᵣ += _fr1_exp(δ,τ,lnδ,lnτ,αᵣ,n_exp,t_exp,d_exp,l)

    #Gaussian-bell-shaped terms
    length(k_gauss) == 0 && return αᵣ
    n_gauss = view(n,k_gauss)
    t_gauss = view(t,k_gauss)
    d_gauss = view(d,k_gauss)
    αᵣ += _fr1_gauss(δ,τ,lnδ,lnτ,αᵣ,n_gauss,t_gauss,d_gauss,η,β,γ,ε)

    #association terms (new)
    length(k_assoc) == 0 && return αᵣ
    lb = length(b)
    lη = length(η)
    k_assoc2 = (lη - lb + 1):lη
    n_assoc = view(n,k_assoc)
    t_assoc = view(t,k_assoc)
    d_assoc = view(d,k_assoc)
    η_assoc = view(η,k_assoc2)
    β_assoc = view(β,k_assoc2)
    γ_assoc = view(γ,k_assoc2)
    ε_assoc = view(ε,k_assoc2)
    b_assoc = b
    αᵣ += _fr1_assoc(δ,τ,lnδ,lnτ,αᵣ,n_assoc,t_assoc,d_assoc,η_assoc,β_assoc,γ_assoc,ε_assoc,b_assoc)
    return αᵣ
end

_frx(model::EmpiricSingleFluid{Nothing},δ,τ) = 0.0

function a_ideal(model::IdealEmpiricSingleFluid,V,T,z=SA[1.])
    Tc = model.properties.Tc
    rhoc = model.properties.rhoc
    N = only(z)
    rho = (N/V)
    δ = rho/rho_c
    τ = Tc/T
    return  _f0(model,δ,τ)
end

a_ideal(model::EmpiricSingleFluid,V,T,z=SA[1.]) = a_ideal(idealmodel(model),V,T,z)

function a_res(model::EmpiricSingleFluid,V,T,z=SA[1.])
    Tc = model.properties.Tc
    rhoc = model.properties.rhoc
    N = only(z)
    rho = (N/V)
    δ = rho/rhoc
    τ = Tc/T
    return  _fr1(model,δ,τ) + _frx(model,δ,τ)
end

function eos(model::EmpiricSingleFluid, V, T, z=SA[1.0])
    R = R_gas(model)
    Tc = model.properties.Tc
    rhoc = model.properties.rhoc
    N = only(z)
    rho = (N/V)
    δ = rho/rhoc
    τ = Tc/T
    return N*R*T*(_f0(model,δ,τ)+_fr1(model,δ,τ) + _frx(model,δ,τ))
end

function eos_res(model::EmpiricSingleFluid,V,T,z=SA[1.0])
    R = R_gas(model)
    Tc = model.consts.Tc
    rhoc = model.consts.rhoc
    N = only(z)
    rho = (N/V)
    δ = rho/rho_c
    τ = Tc/T
    return N*R*T*_fr1(model,δ,τ) + _frx(model,δ,τ)
end

mw(model::EmpiricSingleFluid) = SA[model.properties.Mw]

molecular_weight(model::EmpiricSingleFluid,z = @SVector [1.]) = model.consts.Mw*0.001

T_scale(model::EmpiricSingleFluid,z=SA[1.0]) = model.properties.Tc

p_scale(model::EmpiricSingleFluid,z=SA[1.0]) = model.properties.Pc

lb_volume(model::EmpiricSingleFluid,z=SA[1.0]) = model.properties.lb_volume #finally, an eos model that mentions it max density.

Base.length(::EmpiricSingleFluid) = 1

function Base.show(io::IO,mime::MIME"text/plain",model::EmpiricSingleFluid)
    print(io,"Reference Equation of state for $(model.components[1])")
end

function x0_sat_pure(model::EmpiricSingleFluid,T,z=SA[1.0])
    vv = volume(model.ancilliaries.gas,0.0,T,z)
    vl = volume(model.ancilliaries.liquid,0.0,T,z)
    return (vl,vv)
end

function x0_volume_liquid(model::EmpiricSingleFluid,T,z = SA[1.0])
    volume(model.ancilliaries.liquid,0.0,min(T,model.consts.T_c*one(T)),z)
end

x0_psat(model::EmpiricSingleFluid,T,crit=nothing) = saturation_pressure(model.ancilliaries.saturation,T,SaturationCorrelation())[1]

function x0_saturation_temperature(model::EmpiricSingleFluid,p,z=SA[1.0])
    T = saturation_temperature(model.ancilliaries.saturation,p,SaturationCorrelation())[1]
    vl,vv = x0_sat_pure(model,T)
    return (vl,vv)
end

function crit_pure(model::EmpiricSingleFluid)
    return (model.consts.Tc,model.consts.Pc,1/model.consts.rhoc)
end