
include("structs.jl")

struct EmpiricSingleFluid{𝔸} <: EmpiricHelmholtzModel
    components::Vector{String}
    properties::ESFProperties
    ancilliaries::𝔸
    ideal::ESFIdealParam
    residual::ESFResidualParam
    references::Vector{String}
end

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
EmpiricSingleFluid

struct IdealEmpiricSingleFluid <: IdealModel
    components::Vector{String}
    properties::EmpiricSingleFluidProperties
    ideal::EmpiricSingleFluidIdealParam
    references::Vector{String}
end

function recombine_impl!(model::EmpiricSingleFluid)
    _calc_iterators!(model.residual)
    return model
end

function IdealEmpiricSingleFluid(model::EmpiricSingleFluid)
    return IdealEmpiricSingleFluid(model.components,model.properties,model.ideal,model.references)
end

idealmodel(model::EmpiricSingleFluid) = IdealEmpiricSingleFluid(model)

R_gas(model::EmpiricSingleFluid) = model.properties.Rgas
R_gas(model::IdealEmpiricSingleFluid) = model.properties.Rgas

function _f0(model::Union{EmpiricSingleFluid,IdealEmpiricSingleFluid},δ,τ)
    a₁ = model.ideal.a1
    a₂ = model.ideal.a2
    c₀ = model.ideal.c0
    logδ = log(δ)
    logτ = log(τ)
    α₀ = logδ + a₁ + a₂*τ + c₀*logτ
    n = model.ideal.n_gpe
    #Generalized Plank-Einstein terms
    
    if length(n) != 0
        t = model.ideal.t_gpe
        c = model.ideal.c_gpe
        d = model.ideal.d_gpe
        α₀ +=_f0_gpe(τ,logτ,α₀,n,t,c,d)
    end

    #Power terms
    np = model.ideal.n_p
    if length(np) != 0
        tp = model.ideal.t_p
        α₀ +=_f0_power(τ,logτ,α₀,np,tp)
    end

    return α₀
end

function _fr1(model::EmpiricSingleFluid,δ,τ)

    αᵣ = zero(δ+τ)
    lnδ = log(δ)
    lnτ = log(τ)

    ℙ = model.residual
    n,t,d = ℙ.n,ℙ.t,ℙ.d
    k_pol,k_exp,k_gauss = model.residual.iterators

    #strategy for storing.
    #n, t, d, gauss values, always require views
    #l, b does not require views. they are used just once.

    #Polynomial terms
    n_pol = view(n,k_pol)
    t_pol = view(t,k_pol)
    d_pol = view(d,k_pol)
    αᵣ += _fr1_pol(δ,τ,lnδ,lnτ,αᵣ,n_pol,t_pol,d_pol)

    #Exponential terms.
    if length(k_exp) != 0
        l = ℙ.l
        n_exp = view(n,k_exp)
        t_exp = view(t,k_exp)
        d_exp = view(d,k_exp)
        αᵣ += _fr1_exp(δ,τ,lnδ,lnτ,αᵣ,n_exp,t_exp,d_exp,l)
    end
    #Gaussian-bell-shaped terms
    η,β,γ,ε = ℙ.eta,ℙ.beta,ℙ.gamma,ℙ.epsilon
    if length(k_gauss) != 0
        n_gauss = view(n,k_gauss)
        t_gauss = view(t,k_gauss)
        d_gauss = view(d,k_gauss)
        αᵣ += _fr1_gauss(δ,τ,lnδ,lnτ,αᵣ,n_gauss,t_gauss,d_gauss,η,β,γ,ε)
    end

    #Especial terms are stored in structs.

    #gaoB terms
    if ℙ.gao_b.active
        terms = ℙ.gao_b
        n_gao = terms.n
        t_gao = terms.t
        d_gao = terms.d
        η_gao = terms.eta
        β_gao = terms.beta
        γ_gao = terms.gamma
        ε_gao = terms.epsilon
        b_gao = terms.b
        αᵣ += _fr1_gao(δ,τ,lnδ,lnτ,αᵣ,n_gao,t_gao,d_gao,η_gao,β_gao,γ_gao,ε_gao,b_gao)
    end

    #Non-analytical terms
    if ℙ.na.active
        NA = ℙ.na
        A,B,C,D,aa,bb,ββ,nn = NA.A,NA.B,NA.C,NA.D,NA.a,NA.b,NA.beta,NA.n
        αᵣ += _fr1_na(δ,τ,lnδ,lnτ,αᵣ,A,B,C,D,aa,bb,ββ,nn)
        #αᵣ += iapws95_f0(δ,τ)
    end
 
    return αᵣ
end

function a_ideal(model::IdealEmpiricSingleFluid,V,T,z=SA[1.])
    Tc = model.properties.Tc
    rhoc = model.properties.rhoc
    N = only(z)
    rho = (N/V)
    δ = rho/rhoc
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
    return _fr1(model,δ,τ)
end

function eos(model::EmpiricSingleFluid, V, T, z=SA[1.0])
    R = R_gas(model)
    Tc = model.properties.Tc
    rhoc = model.properties.rhoc
    N = only(z)
    rho = (N/V)
    δ = rho/rhoc
    τ = Tc/T
    return N*R*T*(_f0(model,δ,τ)+_fr1(model,δ,τ))
end

function eos_res(model::EmpiricSingleFluid,V,T,z=SA[1.0])
    R = R_gas(model)
    Tc = model.properties.Tc
    rhoc = model.properties.rhoc
    N = only(z)
    rho = (N/V)
    δ = rho/rho_c
    τ = Tc/T
    return N*R*T*_fr1(model,δ,τ)
end

mw(model::EmpiricSingleFluid) = SA[model.properties.Mw]

molecular_weight(model::EmpiricSingleFluid,z = @SVector [1.]) = model.properties.Mw*0.001

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
    volume(model.ancilliaries.liquid,0.0,min(T,model.properties.T_c*one(T)),z)
end

x0_psat(model::EmpiricSingleFluid,T,crit=nothing) = saturation_pressure(model.ancilliaries.saturation,T,SaturationCorrelation())[1]

function x0_saturation_temperature(model::EmpiricSingleFluid,p,z=SA[1.0])
    T = saturation_temperature(model.ancilliaries.saturation,p,SaturationCorrelation())[1]
    vl,vv = x0_sat_pure(model,T)
    return (vl,vv)
end

function crit_pure(model::EmpiricSingleFluid)
    Tc = model.properties.Tc
    Vc = 1/model.properties.rhoc
    Pc = model.properties.Pc

    return (Tc,Pc,Vc)
end

function tryparse_units(val,unit)
    result = try
        unit_parsed = Unitful.uparse(unit)
        ThermoState.normalize_units(val*unit_parsed)
    catch
        val
    end
    return result
end

function fff(path::String)
    _path = only(flattenfilepaths(String[],path))

    json_string = read(_path, String)
    data = JSON3.read(json_string)
end

include("parser.jl")