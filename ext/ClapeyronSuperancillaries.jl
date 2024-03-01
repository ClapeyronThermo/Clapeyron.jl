module ClapeyronSuperancillaries
using EoSSuperancillaries
using Clapeyron
const C = Clapeyron

#(:PCSAFT,:PPCSAFT,:QPPCSAFT),
const SuperancPCSAFT = Union{C.PCSAFT,C.PPCSAFT,C.QPPCSAFT,C.pharmaPCSAFT}

function can_superanc(model::SuperancPCSAFT)
    val = true
    val = val & (length(model.params.epsilon_assoc.values.values) == 0)
    if hasfield(typeof(model.params),:dipole)
        val = val & all(iszero,model.params.dipole.values)
    end
    return val && C.SUPERANC_ENABLED[]
end

x0_sat_pure_default(model::SuperancPCSAFT,T) = Clapeyron.x0_sat_pure_virial(model,T)

function Δσ(model,T)
    if model isa C.pharmaPCSAFT
        k = C.water08_k(model)
        if k > 0 
            σ += C.Δσh20(T)
        end
    else
        return zero(T)
    end
end

function get_pcsaft_consts(model)
    m = model.params.segment.values[1]
    ϵ = model.params.epsilon.values[1]
    σ = model.params.sigma.values[1]
    return m,ϵ,σ
end

function C.x0_sat_pure(model::SuperancPCSAFT,T)
    can_superanc(model) || return x0_sat_pure_default(model,T)
    m,ϵ,σ = get_pcsaft_consts(model)
    σ += Δσ(model,T)
    T̃ = T/ϵ
    θ,status = pcsaft_theta(T̃,m)
    if status == :over_Tmax
        return θ,θ
    elseif status != :inrange
        return x0_sat_pure_default(model,T)
    else
        ρ̃l,ρ̃v = pcsaft_rhosat_reduced(θ,m)
        N_Aσ3 = C.N_A*σ*σ*σ
        return N_Aσ3/ρ̃l,N_Aσ3/ρ̃v
    end
end

function x0_crit_pure_default(model::SuperancPCSAFT) 
    lb_v = C.lb_volume(model)
    (2.0, log10(lb_v/0.3))
end

function C.x0_crit_pure(model::SuperancPCSAFT)
    can_superanc(model) || return x0_crit_pure_default(model)
    m,ϵ,σ = get_pcsaft_consts(model)
    if 1.0 <= m <= 64.0
        Tc = pcsaft_tc(m,ϵ)
        vc = pcsaft_vc(m,σ + Δσ(model,Tc))
        return Tc/ϵ,log10(vc)
    else
        return x0_crit_pure_default(model)
    end
end

function C.crit_pure(model::SuperancPCSAFT)
    can_superanc(model) || return C.crit_pure(model,C.x0_crit_pure(model))
    m,ϵ,σ = get_pcsaft_consts(model)
    if 1.0 <= m <= 64.0
        Tc = pcsaft_tc(m,ϵ)
        vc = pcsaft_vc(m,σ + Δσ(model,Tc))
        pc = pressure(model,vc,Tc)
        return Tc,pc,vc
    else
        return crit_pure(model,C.x0_crit_pure_default(model))
    end
end

function C.saturation_pressure_impl(model::SuperancPCSAFT,T,method::C.SuperAncSaturation)
    _0 = zero(T + 1.0 + oneunit(eltype(model)))
    nan = _0/_0
    fail = (nan,nan,nan)
    can_superanc(model) || return fail
    m,ϵ,σ = get_pcsaft_consts(model)
    σ += Δσ(model,T)
    vl,vv = pcsaft_vsat(T,m,ϵ,σ)
    p = pressure(model,vv,T)
    return p,vl,vv
end

const SuperancCubic = Union{C.vdW,C.PR,C.RK}

x0_sat_pure_default(model::SuperancCubic,T) = Clapeyron.x0_sat_pure_cubic_ab(model,T)

function can_superanc(model::SuperancCubic)
    return C.SUPERANC_ENABLED[]
end

function C.x0_sat_pure(model::SuperancCubic,T)
    can_superanc(model) || return x0_sat_pure_default(model,T)
    a,b,c = C.cubic_ab(model,1e-3,T,C.SA[1.0])
    ac,bc = model.params.a[1],model.params.b[1]
    _0 = zero(a)
    Tc = model.params.Tc.values[1]
    if T > Tc
        nan = _0/_0
        return nan,nan
    end
    T̃ = T*Rgas(model)*b/a
    T̃c = Tc*Rgas(model)*bc/ac
    T̃ < 0.1*T̃c && return x0_sat_pure_default(model,T)
    if model isa C.vdW
        return vdw_vsat(T,a,b) .- c
    elseif model isa C.RK
        return rk_vsat(T,a,b) .- c
    else #model isa C.PR
        return pr_vsat(T,a,b) .- c
    end
end

function __init__()
    C.use_superancillaries!(true)
end

end #module
