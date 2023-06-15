
include("structs.jl")

#term dispatch. function definitions are in term_functions.jl

function a_term(term::NonAnalyticTerm,δ,τ,lnδ,lnτ,_0)
    if term.active    
        A,B,C,D,a,b,β,n = term.A,term.B,term.C,term.D,term.a,term.b,term.beta,term.n
        αᵣ = term_ar_na(δ,τ,lnδ,lnτ,_0,A,B,C,D,a,b,β,n)
    else
        αᵣ = _0
    end
    return αᵣ
end

function a_term(term::GaoBTerm,δ,τ,lnδ,lnτ,_0)
    if term.active    
        n = term.n
        t = term.t
        d = term.d
        η = term.eta
        β = term.beta
        γ = term.gamma
        ε = term.epsilon
        b = term.b
        αᵣ = term_ar_gaob(δ,τ,lnδ,lnτ,_0,n,t,d,η,β,γ,ε,b)
    else
        αᵣ = _0
    end
    return αᵣ
end

function a_term(term::Associating2BTerm,δ,τ,lnδ,lnτ,_0)
    if term.active    
        ε = term.epsilonbar
        κ = term.kappabar
        a = term.a
        m = term.m
        v̄ₙ = term.vbarn
        αᵣ = term_ar_assoc2b(δ,τ,lnδ,lnτ,_0,ε,κ,a,m,v̄ₙ)
    else
        αᵣ = _0
    end
    return αᵣ
end

struct EmpiricSingleFluid{𝔸} <: EmpiricHelmholtzModel
    components::Vector{String}
    properties::ESFProperties
    ancillaries::𝔸
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

reduced_a_ideal(model::EmpiricSingleFluid,τ) = reduced_a_ideal(model.ideal,τ)
reduced_a_ideal(model::IdealEmpiricSingleFluid,τ) = reduced_a_ideal(model.ideal,τ)

function reduced_a_ideal(model::EmpiricSingleFluidIdealParam,τ)
    a₁ = model.a1
    a₂ = model.a2
    c₀ = model.c0
    logτ = log(τ)
    α₀ = a₁ + a₂*τ + c₀*logτ
   
    #Generalized Plank-Einstein terms
    n = model.n_gpe
    if length(n) != 0
        t = model.t_gpe
        c = model.c_gpe
        d = model.d_gpe
        α₀ +=term_a0_gpe(τ,logτ,α₀,n,t,c,d)
    end

    #Power terms
    np = model.n_p
    if length(np) != 0
        tp = model.ideal.t_p
        α₀ += term_a0_power(τ,logτ,α₀,np,tp)
    end

    #GERG-2008 terms
    n_gerg = model.n_gerg
    if length(n_gerg) != 0
        v_gerg = model.v_gerg
        α₀ += term_a0_gerg2008(τ,logτ,α₀,n_gerg,v_gerg)
    end

    return α₀
end

reduced_a_res(model::EmpiricSingleFluid,δ,τ,lnδ = log(δ),lnτ = log(τ)) = reduced_a_res(model.residual,δ,τ,lnδ,lnτ)

function reduced_a_res(model::EmpiricSingleFluidResidualParam,δ,τ,lnδ = log(δ),lnτ = log(τ))
    _0 = zero(δ+τ)
    αᵣ = _0
    ℙ = model
    n,t,d = ℙ.n,ℙ.t,ℙ.d
    k_pol,k_exp,k_gauss = ℙ.iterators

    #strategy for storing.
    #n, t, d, gauss values, always require views
    #l, b does not require views. they are used just once.

    #Polynomial terms
    n_pol = view(n,k_pol)
    t_pol = view(t,k_pol)
    d_pol = view(d,k_pol)
    αᵣ += term_ar_pol(δ,τ,lnδ,lnτ,αᵣ,n_pol,t_pol,d_pol)

    #Exponential terms.
    if length(k_exp) != 0
        l,g = ℙ.l,ℙ.g
        n_exp = view(n,k_exp)
        t_exp = view(t,k_exp)
        d_exp = view(d,k_exp)
        αᵣ += term_ar_exp(δ,τ,lnδ,lnτ,αᵣ,n_exp,t_exp,d_exp,l,g)
    end

    #Gaussian bell-shaped terms
    η,β,γ,ε = ℙ.eta,ℙ.beta,ℙ.gamma,ℙ.epsilon
    if length(k_gauss) != 0
        n_gauss = view(n,k_gauss)
        t_gauss = view(t,k_gauss)
        d_gauss = view(d,k_gauss)
        αᵣ += term_ar_gauss(δ,τ,lnδ,lnτ,αᵣ,n_gauss,t_gauss,d_gauss,η,β,γ,ε)
    end

    #Especial terms are stored in structs.

    #gaoB terms
    αᵣ += a_term(ℙ.gao_b,δ,τ,lnδ,lnτ,_0)
    
    #Non-analytical terms
    αᵣ += a_term(ℙ.na,δ,τ,lnδ,lnτ,_0)

    #associating terms.
    αᵣ += a_term(ℙ.assoc,δ,τ,lnδ,lnτ,_0)
 
    return αᵣ
end

function __get_k_alpha0(model)
    R0 = model.ideal.R0
    if iszero(R0)
        return 1.0
    else
        R = model.properties.Rgas
        return R0/R
    end
end

function a_ideal(model::IdealEmpiricSingleFluid,V,T,z=SA[1.],k = __get_k_alpha0(model))
    Tc = model.properties.Tc
    rhoc = model.properties.rhoc
    N = only(z)
    rho = (N/V)
    δ = rho/rhoc
    τ = Tc/T
    α0 = reduced_a_ideal(model,τ)
    return k*α0 + log(δ)
end

a_ideal(model::EmpiricSingleFluid,V,T,z=SA[1.]) = a_ideal(idealmodel(model),V,T,z)

function a_res(model::EmpiricSingleFluid,V,T,z=SA[1.])
    Tc = model.properties.Tc
    rhoc = model.properties.rhoc
    N = only(z)
    rho = (N/V)
    δ = rho/rhoc
    τ = Tc/T
    return reduced_a_res(model,δ,τ)
end

function eos(model::EmpiricSingleFluid, V, T, z=SA[1.0])
    R = R_gas(model)
    Tc = model.properties.Tc
    rhoc = model.properties.rhoc
    N = only(z)
    rho = (N/V)
    δ = rho/rhoc
    τ = Tc/T
    k = __get_k_alpha0(model)
    return N*R*T*(log(δ) + k*reduced_a_ideal(model,τ) + reduced_a_res(model,δ,τ))
end

function eos_res(model::EmpiricSingleFluid,V,T,z=SA[1.0])
    R = R_gas(model)
    Tc = model.properties.Tc
    rhoc = model.properties.rhoc
    N = only(z)
    rho = (N/V)
    δ = rho/rhoc
    τ = Tc/T
    return N*R*T*reduced_a_res(model,δ,τ)
end

mw(model::EmpiricSingleFluid) = SA[model.properties.Mw]

molecular_weight(model::EmpiricSingleFluid,z = @SVector [1.]) = model.properties.Mw*0.001

T_scale(model::EmpiricSingleFluid,z=SA[1.0]) = model.properties.Tc

p_scale(model::EmpiricSingleFluid,z=SA[1.0]) = model.properties.Pc

lb_volume(model::EmpiricSingleFluid,z=SA[1.0]) = model.properties.lb_volume #finally, an eos model that mentions it max density.

Base.length(::EmpiricSingleFluid) = 1

function Base.show(io::IO,mime::MIME"text/plain",model::EmpiricSingleFluid)
    println(io,"MultiParameter Equation of state for $(model.components[1]):")
    show_multiparameter_coeffs(io,model.residual)
end

function Base.show(io::IO,mime::MIME"text/plain",model::IdealEmpiricSingleFluid)
    println(io,"Ideal MultiParameter Equation of state for $(model.components[1]):")
    show_multiparameter_coeffs(io,model.ideal)

end

function x0_sat_pure(model::EmpiricSingleFluid,T,z=SA[1.0])
    vv = volume(model.ancillaries.gas,0.0,T,z)
    vl = x0_volume_liquid(model,T,z)
    return (vl,vv)
end

function x0_volume_liquid(model::EmpiricSingleFluid,T,z = SA[1.0])
    lb_v = lb_volume(model)
    vl_tp = 1/model.properties.rhol_tp
    vl_anc = volume(model.ancillaries.liquid,0.0,min(T,model.properties.Tc*one(T)),z)
    isnan(vl_tp) && (vl_tp = 0.0)
    isnan(vl_anc) && (vl_anc = 0.0)
    return max(vl_tp,vl_anc,1.01*lb_v)
end

x0_psat(model::EmpiricSingleFluid,T,crit=nothing) = saturation_pressure(model.ancillaries.saturation,T,SaturationCorrelation())[1]

function x0_saturation_temperature(model::EmpiricSingleFluid,p)
    T = saturation_temperature(model.ancillaries.saturation,p,SaturationCorrelation())[1]
    vl,vv = x0_sat_pure(model,T)
    return (T,vl,vv)
end

function crit_pure(model::EmpiricSingleFluid)
    Tc = model.properties.Tc
    Vc = 1/model.properties.rhoc
    Pc = model.properties.Pc
    return (Tc,Pc,Vc)
end

include("parser.jl")