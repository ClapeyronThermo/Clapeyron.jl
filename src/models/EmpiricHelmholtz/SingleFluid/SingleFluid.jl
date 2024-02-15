
include("structs.jl")

const EmpiricAncillary = CompositeModel{FluidCorrelation{PolExpVapour, PolExpLiquid, PolExpSat}, Nothing}
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

struct SingleFluid{𝔸} <: EmpiricHelmholtzModel
    components::Vector{String}
    properties::ESFProperties
    ancillaries::𝔸
    ideal::ESFIdealParam
    residual::ESFResidualParam
    references::Vector{String}
end

struct SingleFluidIdeal <: IdealModel
    components::Vector{String}
    properties::SingleFluidProperties
    ideal::SingleFluidIdealParam
    references::Vector{String}
end

function recombine_impl!(model::SingleFluid)
    _calc_iterators!(model.residual)
    return model
end

function SingleFluidIdeal(model::SingleFluid)
    return SingleFluidIdeal(model.components,model.properties,model.ideal,model.references)
end

idealmodel(model::SingleFluid) = SingleFluidIdeal(model)

R_gas(model::SingleFluid) = model.properties.Rgas
R_gas(model::SingleFluidIdeal) = model.properties.Rgas

reduced_a_ideal(model::SingleFluid,τ) = reduced_a_ideal(model.ideal,τ)
reduced_a_ideal(model::SingleFluidIdeal,τ) = reduced_a_ideal(model.ideal,τ)

function reduced_a_ideal(model::SingleFluidIdealParam,τ)
    a₁ = model.a1
    a₂ = model.a2
    c₀ = model.c0
    c₁ = model.c1
    logτ = log(τ)
    α₀ = a₁ + a₂*τ + c₀*logτ
    if !iszero(c₁)
        α₀ += c₁*τ*logτ
    end
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
        tp = model.t_p
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

reduced_a_res(model::SingleFluid,δ,τ,lnδ = log(δ),lnτ = log(τ)) = reduced_a_res(model.residual,δ,τ,lnδ,lnτ)

function reduced_a_res(ℙ::MultiParameterParam,δ,τ,lnδ = log(δ),lnτ = log(τ))
    _0 = zero(δ+τ)
    αᵣ = _0
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
    __has_extra_params(ℙ) || return αᵣ

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

function __set_Rgas(pure,Rgas)
    components,p,ancillaries,ideal,residual,references = pure.components,pure.properties,pure.ancillaries,pure.ideal,pure.residual,pure.references
    Mw,Tr,rhor,lb_volume,Tc,Pc,rhoc,Ttp,ptp,rhov_tp,rhol_tp,acentricfactor = p.Mw, p.Tr, p.rhor, p.lb_volume, p.Tc, p.Pc, p.rhoc, p.Ttp, p.ptp, p.rhov_tp, p.rhol_tp, p.acentricfactor
    properties = ESFProperties(Mw,Tr,rhor,lb_volume,Tc,Pc,rhoc,Ttp,ptp,rhov_tp,rhol_tp,acentricfactor,Rgas)
    return SingleFluid(components,properties,ancillaries,ideal,residual,references)
end

function a_ideal(model::SingleFluidIdeal,V,T,z=SA[1.],k = __get_k_alpha0(model))
    Tc = model.properties.Tc
    rhoc = model.properties.rhoc
    N = sum(z)
    τ = Tc/T
    α0 = reduced_a_ideal(model,τ)
    #this form separates the dependency of z from the dependency of V, allowing for precise derivatives of the ideal part.
    logδ = log(N/rhoc) - log(V)
    return k*α0 + logδ
end

v_scale(model::SingleFluid,z = SA[1.0],∑z = sum(z)) = 1/∑z/model.properties.rhoc
v_scale(model::SingleFluidIdeal,z = SA[1.0],∑z = sum(z)) = 1/∑z/model.properties.rhoc

a_ideal(model::SingleFluid,V,T,z=SA[1.]) = a_ideal(idealmodel(model),V,T,z)

function a_res(model::SingleFluid,V,T,z=SA[1.])
    Tc = model.properties.Tc
    rhoc = model.properties.rhoc
    N = sum(z)
    δ = N/(rhoc*V)
    τ = Tc/T
    return reduced_a_res(model,δ,τ)
end

function eos(model::SingleFluid, V, T, z=SA[1.0])
    R = R_gas(model)
    Tc = model.properties.Tc
    rhoc = model.properties.rhoc
    N = sum(z)
    δ = N/(rhoc*V)
    τ = Tc/T
    k = __get_k_alpha0(model)
    logδ = log(δ)
    ref_a = model.ideal.ref_a
    a0,a1 = ref_a[1],ref_a[2]
    return N*R*T*(logδ + k*reduced_a_ideal(model,τ) + reduced_a_res(model,δ,τ)) + N*(a0 + a1*T)
end

function eos_res(model::SingleFluid,V,T,z=SA[1.0])
    R = R_gas(model)
    Tc = model.properties.Tc
    rhoc = model.properties.rhoc
    N = sum(z)
    δ = N/(rhoc*V)
    τ = Tc/T
    return N*R*T*reduced_a_res(model,δ,τ)
end

mw(model::SingleFluid) = SA[model.properties.Mw]

molecular_weight(model::SingleFluid,z = @SVector [1.]) = model.properties.Mw*0.001

T_scale(model::SingleFluid,z=SA[1.0]) = model.properties.Tc

p_scale(model::SingleFluid,z=SA[1.0]) = model.properties.Pc

lb_volume(model::SingleFluid,z=SA[1.0]) = model.properties.lb_volume #finally, an eos model that mentions it max density.

Base.length(::SingleFluid) = 1

function Base.show(io::IO,mime::MIME"text/plain",model::SingleFluid)
    println(io,"MultiParameter Equation of state for $(model.components[1]):")
    show_multiparameter_coeffs(io,model.residual)
end

function Base.show(io::IO,mime::MIME"text/plain",model::SingleFluidIdeal)
    println(io,"Ideal MultiParameter Equation of state for $(model.components[1]):")
    show_multiparameter_coeffs(io,model.ideal)
end

function x0_sat_pure(model::SingleFluid,T)
    z=SA[1.0]
    gas_ancillary = model.ancillaries.fluid.gas
    vv = volume(gas_ancillary,0.0,T,z)
    vl = x0_volume_liquid(model,T,)
    return (vl,vv)
end

function x0_volume_liquid(model::SingleFluid,T,z = SA[1.0])
    lb_v = lb_volume(model)
    vl_tp = 1/model.properties.rhol_tp
    liquid_ancillary = model.ancillaries.fluid.liquid
    vl_anc = volume(liquid_ancillary,0.0,min(T,model.properties.Tc*one(T)),z)
    isnan(vl_tp) && (vl_tp = 0.0)
    isnan(vl_anc) && (vl_anc = 0.0)
    return max(vl_tp,vl_anc,1.01*lb_v)
end

x0_psat(model::SingleFluid,T,crit=nothing) = saturation_pressure(model.ancillaries.fluid.saturation,T,SaturationCorrelation())[1]

function x0_saturation_temperature(model::SingleFluid,p)
    T = saturation_temperature(model.ancillaries.fluid.saturation,p,SaturationCorrelation())[1]
    vl,vv = x0_sat_pure(model,T)
    return (T,vl,vv)
end

function crit_pure(model::SingleFluid)
    Tc = model.properties.Tc
    Vc = 1/model.properties.rhoc
    Pc = model.properties.Pc
    return (Tc,Pc,Vc)
end

include("parser.jl")

export EmpiricAncillary, SingleFluid, SingleFluidIdeal
