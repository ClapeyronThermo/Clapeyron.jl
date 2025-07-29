
abstract type MultiParameterTerm end

include("terms/polexpgauss.jl")
include("terms/exp2.jl")
include("terms/gaob.jl")
include("terms/nonanalytic.jl")
include("terms/assoc2B.jl")

include("structs.jl")

const EmpiricAncillary = CompositeModel{FluidCorrelation{PolExpVapour, PolExpLiquid, PolExpSat, Nothing}, Nothing}
#term dispatch. function definitions are in term_functions.jl

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

is_splittable(model::SingleFluid) = false
is_splittable(model::SingleFluidIdeal) = false

function recombine_impl!(model::SingleFluid)
    _calc_iterators!(model.residual)
    return model
end

function SingleFluidIdeal(model::SingleFluid)
    return SingleFluidIdeal(model.components,model.properties,model.ideal,model.references)
end

idealmodel(model::SingleFluid) = SingleFluidIdeal(model)

Rgas(model::SingleFluid) = model.properties.Rgas
Rgas(model::SingleFluidIdeal) = model.properties.Rgas

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

    #pol+exp+gauss terms
    αᵣ += a_term(ℙ.polexpgauss,δ,τ,lnδ,lnτ,_0)

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

function a_ideal(model::SingleFluidIdeal,V,T,z=SA[1.],k = __get_k_alpha0(model))
    Tc = model.properties.Tr
    rhoc = model.properties.rhor
    N = sum(z)
    τ = Tc/T
    α0 = reduced_a_ideal(model,τ)
    #this form separates the dependency of z from the dependency of V, allowing for precise derivatives of the ideal part.
    logδ = log(N/rhoc) - log(V)
    return k*α0 + logδ
end

v_scale(model::SingleFluid,z = SA[1.0],∑z = sum(z)) = 1/∑z/model.properties.rhor
v_scale(model::SingleFluidIdeal,z = SA[1.0],∑z = sum(z)) = 1/∑z/model.properties.rhor

a_ideal(model::SingleFluid,V,T,z=SA[1.]) = a_ideal(idealmodel(model),V,T,z)

function a_res(model::SingleFluid,V,T,z=SA[1.])
    Tc = model.properties.Tr
    rhoc = model.properties.rhor
    N = sum(z)
    δ = N/(rhoc*V)
    τ = Tc/T
    return reduced_a_res(model,δ,τ)
end

function eos_impl(model::SingleFluid, V, T, z)
    R = Rgas(model)
    Tc = model.properties.Tr
    rhoc = model.properties.rhor
    N = sum(z)
    δ = N/(rhoc*V)
    τ = Tc/T
    k = __get_k_alpha0(model)
    logδ = log(δ)
    ref_a = model.ideal.ref_a
    a0,a1 = ref_a[1],ref_a[2] #reference state evaluation
    return N*R*T*(logδ + k*reduced_a_ideal(model,τ) + reduced_a_res(model,δ,τ,logδ)) + N*(a0 + a1*T)
end

mw(model::SingleFluid) = SA[model.properties.Mw]

T_scale(model::SingleFluid,z) = model.properties.Tc

p_scale(model::SingleFluid,z) = model.properties.Pc

lb_volume(model::SingleFluid,z) = model.properties.lb_volume #finally, an eos model that mentions it max density.

Base.length(::SingleFluid) = 1

function Base.show(io::IO,mime::MIME"text/plain",model::SingleFluid)
    println(io,"MultiParameter Equation of state for $(model.components[1]):")
    show_multiparameter_coeffs(io,model.residual)
end

function Base.show(io::IO,mime::MIME"text/plain",model::SingleFluidIdeal)
    println(io,"Ideal MultiParameter Equation of state for $(model.components[1]):")
    show_multiparameter_coeffs(io,model.ideal)
end
has_fast_crit_pure(model::SingleFluid) = true

function x0_sat_pure(model::SingleFluid,T)
    z=SA[1.0]
    Ttp0 = model.properties.Ttp*one(T)
    gas_ancillary = model.ancillaries.fluid.gas
    liquid_ancillary = model.ancillaries.fluid.liquid
    if isfinite(Ttp0) & (Ttp0 > 0) #we have triple point information:
        Ttp = Ttp0
    else #we suppose triple point = 0.4 Tc
        Ttp = 0.4*model.properties.Tc*one(T)
    end
    if T > model.properties.Tc
        _0 = zero(Base.promote_eltype(model,T))
        _nan = _0/_0
        return _nan,_nan
    elseif T >= Ttp
        vv = volume(gas_ancillary,0.0,T,z)
        vl = volume(liquid_ancillary,0.0,T,z)
        return vl,vv
    else #we know that T < Ttp
        vvtp = volume(gas_ancillary,0.0,Ttp,z)
        vltp = volume(liquid_ancillary,0.0,Ttp,z)
        ptp = pressure(model,vvtp,Ttp)
        R = Rgas(model)
        dpdT = dpdT_saturation(model,vvtp,vltp,Ttp)
        dTinvdlnp = -ptp/(dpdT*T*T)
        ΔTinv = 1/T - 1/Ttp
        psat = exp(ΔTinv/dTinvdlnp)*ptp
        vv = Rgas(model)*T/psat
        vl = x0_volume_liquid(model,psat,T)
        return vl,vv
    end
end

function x0_volume_liquid_lowT(model::SingleFluid,p,T,z)
    _1 = one(Base.promote_eltype(model,p,T,z))
    lb_v = lb_volume(model,T,z)*_1
    vl_lbv = 1.01*lb_v
    ancillary = model.ancillaries
    Tc = model.properties.Tc
    Ttp = model.properties.Ttp

    if Ttp < T < Tc
        vᵢ = volume(ancillary,p,T,z,phase = :l)
        return vᵢ
    elseif Ttp < T
        vᵢ = volume(ancillary,p,Ttp,z,phase = :l)
        pvi = pressure(model,vᵢ,T)
        pp = max(_1*p,pvi)
        Tᵢ = _1*Ttp
        #chill from p,Ttp to p,T
        return volume_chill(model,pp,T,z,vᵢ,Tᵢ)
    else #T > Tc and p < pc, this is gas-like supercritical fluid
        #this is always a gas volume, so starting from the lowest volume does not hurt
        return vl_lbv
    end
end

function x0_volume_liquid(model::SingleFluid,p,T,z)
    _1 = one(Base.promote_eltype(model,p,T,z))
    lb_v = lb_volume(model,T,z)*_1
    vl_lbv = 1.001*lb_v
    Tc = model.properties.Tc
    Pc = model.properties.Pc
    Ttp = model.properties.Ttp
    ptp = model.properties.ptp
    (!isfinite(Ttp) | (Ttp < 0)) && (Ttp = 0.4*Tc)
    (!isfinite(ptp) | (ptp < 0)) && (ptp = zero(ptp))
    if p > Pc
        #supercritical conditions, liquid
        #https://doi.org/10.1016/j.ces.2018.08.043 gives an aproximation of the pv curve at T = Tc
        #=
        abs(1 - P/Pc) = abs(1-Vc/V)^(1/Zc)
        abs(1 - P/Pc)^Zc = abs(1-Vc/V)
        if T > Tc then V > V(T = Tc) ≈ V_crit(P)
        =#
        vc = 1/model.properties.rhoc
        pc = model.properties.Pc
        Zc = pc*vc/(Rgas(model)*Tc)
        ΔVrm1 = _1*(abs(1 - p/pc))^Zc
        v_crit_aprox = vc/(ΔVrm1 + 1)
        vhi =  max(vl_lbv,v_crit_aprox)
        phi = pressure(model,vhi,T,z)
        #we suppose that V < Vc (liquid state), then the volume solver converges really well with this initial guess
        if T >= Tc
            if phi > p
                return vhi
            elseif phi <= p <= pressure(model,lb_v,T,z)
                return volume_bracket_refine(model,p,T,z,lb_v,vhi)
            else
                return vl_lbv
            end
        else
            #we want two points: psat-vsat and phi-vhi
            #we can interpolate those to calculate an initial volume
            vsat = x0_volume_liquid_lowT(model,p,T,z)
            psat = pressure(model,vsat,T,z)

            if vhi > vsat
                vhi = 0.9*vsat + 0.1*lb_v
                phi = pressure(model,vhi,T,z)
            end

            if phi <= p
                return volume_bracket_refine(model,p,T,z,vhi,lb_v)
            elseif psat < p < phi
                return volume_bracket_refine(model,p,T,z,vhi,vsat)
            else
                return vsat
            end
        end
    else
        return x0_volume_liquid_lowT(model,p,T,z)
    end
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
