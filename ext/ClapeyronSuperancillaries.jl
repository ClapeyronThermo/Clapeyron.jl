module ClapeyronSuperancillaries
import EoSSuperancillaries
import Clapeyron
const C = Clapeyron
const ES = EoSSuperancillaries
#(:PCSAFT,:PCPSAFT,:QPCPSAFT),
const SuperancPCSAFT = Union{C.PCSAFT,C.PCPSAFT,C.QPCPSAFT,C.pharmaPCSAFT}

function can_superanc(model::SuperancPCSAFT)
    val = true
    val = val & (length(model.params.epsilon_assoc.values.values) == 0)
    if hasfield(typeof(model.params),:dipole)
        val = val & all(iszero,model.params.dipole.values)
    end
    return val && C.SUPERANC_ENABLED[]
end

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

function C.x0_sat_pure(model::SuperancPCSAFT,T,crit = nothing)
    if !can_superanc(model)
        if crit === nothing
            _,vl,vv = C.x0_sat_pure_virial(model,T)
            return vl,vv
        else
            _,vl,vv = C.x0_sat_pure_crit(model,T,crit)
            return vl,vv
        end
    end
    m,ϵ,σ = get_pcsaft_consts(model)
    σ += Δσ(model,T)
    T̃ = T/ϵ
    θ,status = ES.pcsaft_theta(T̃,m)
    if status == :over_Tmax
        return θ,θ
    elseif status == :nonfinite
        nan = zero(θ)/zero(θ)
        return nan,nan
    elseif status == :below_Tmin
        #WARNING: this is really wrong in Float64 arithmetic, but helps in extended precision
        T0 = ES.pcsaft_tc(m,ϵ)*exp(-2.20078778)*m^0.37627892
        ρ̃l,ρ̃v = ES.pcsaft_rhosat_reduced(zero(T),m)
        N_Aσ3 = C.N_A*σ*σ*σ
        vli,vvi = N_Aσ3/ρ̃l,N_Aσ3/ρ̃v
        pii = C.pressure(model,vvi,T0)
        dpdT = C.dpdT_pure(model,vvi,vli,T0)
        dTinvdlnp = -pii/(dpdT*T*T)
        Δlnp = (1/T - 1/T0)/dTinvdlnp
        p = exp(Δlnp)*pii
        return vli, C.Rgas(model)*T/p

    elseif status != :inrange #segment diameter over 64, weird but we should still offer an initial point.
        if crit === nothing
            _,vl,vv = C.x0_sat_pure_virial(model,T)
            return vl,vv
        else
            _,vl,vv = C.x0_sat_pure_crit(model,T,crit)
            return vl,vv
        end
    else
        ρ̃l,ρ̃v = ES.pcsaft_rhosat_reduced(θ,m)
        N_Aσ3 = C.N_A*σ*σ*σ
        return N_Aσ3/ρ̃l,N_Aσ3/ρ̃v
    end
end

function C.x0_crit_pure(model::SuperancPCSAFT)
    if !can_superanc(model) 
        lb_v = C.lb_volume(model)
        return (2.0*oneunit(lb_v), log10(lb_v/0.3))
    end
    m,ϵ,σ = get_pcsaft_consts(model)
    if 1.0 <= m <= 64.0
        Tc = ES.pcsaft_tc(m,ϵ)
        vc = ES.pcsaft_vc(m,σ + Δσ(model,Tc))
        return Tc/ϵ,log10(vc)
    else
        lb_v = C.lb_volume(model)
        return (2.0*oneunit(lb_v), log10(lb_v/0.3))
    end
end

function C.crit_pure(model::SuperancPCSAFT)
    can_superanc(model) || return C.crit_pure(model,C.x0_crit_pure(model))
    m,ϵ,σ = get_pcsaft_consts(model)
    if 1.0 <= m <= 64.0 && eltype(model) === Float64
        Tc = ES.pcsaft_tc(m,ϵ)
        vc = ES.pcsaft_vc(m,σ + Δσ(model,Tc))
        pc = C.pressure(model,vc,Tc)
        return Tc,pc,vc
    else
        return C.crit_pure(model,C.x0_crit_pure(model))
    end
end

function C.x0_psat(model::SuperancPCSAFT,T,crit = nothing)
    can_superanc(model) || return C.x0_psat(model, T, nothing)
    m,ϵ,σ = get_pcsaft_consts(model)
    if 1.0 <= m <= 64.0
        vl = ES.pcsaft_vlsat(T,m,ϵ,σ)
        return C.pressure(model,vl,T)
    else
        return C.x0_psat(model, T, crit)
    end
end

function C.saturation_pressure_impl(model::SuperancPCSAFT,T,method::C.SuperAncSaturation)
    _0 = zero(T + 1.0 + oneunit(eltype(model)))
    nan = _0/_0
    fail = (nan,nan,nan)
    can_superanc(model) || return fail
    m,ϵ,σ = get_pcsaft_consts(model)
    σ += Δσ(model,T)
    vl,vv = ES.pcsaft_vsat(T,m,ϵ,σ)
    p = C.pressure(model,vv,T)
    return p,vl,vv
end

const SuperancCubic = Union{C.vdW,C.PR,C.RK}

function can_superanc(model::SuperancCubic)
    return C.SUPERANC_ENABLED[]
end

function C.x0_sat_pure(model::SuperancCubic,T,crit = nothing)
    can_superanc(model) || return C.x0_sat_pure_cubic_ab(model,T)
    a,b,c = C.cubic_ab(model,1e-3,T,C.SA[1.0])
    ac,bc = model.params.a[1],model.params.b[1]
    _0 = zero(a)
    Tc = model.params.Tc.values[1]
    if T > Tc
        nan = _0/_0
        return nan,nan
    end
    T̃ = T*C.Rgas(model)*b/a
    T̃c = Tc*C.Rgas(model)*bc/ac
    T̃ < 0.1*T̃c && return C.x0_sat_pure_cubic_ab(model,T)
    if model isa C.vdW
        return ES.vdw_vsat(T,a,b) .- c
    elseif model isa C.RK
        return ES.rk_vsat(T,a,b) .- c
    else #model isa C.PR
        return ES.pr_vsat(T,a,b) .- c
    end
end

function __init__()
    C.use_superancillaries!(true)
end

end #module
