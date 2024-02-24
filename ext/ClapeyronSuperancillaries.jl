module ClapeyronSuperancillaries
using EoSSuperancillaries
using Clapeyron
const C = Clapeyron

#(:PCSAFT,:PPCSAFT,:QPPCSAFT),
const SuperancPCSAFT = Union{C.PCSAFT,C.PPCSAFT,C.QPPCSAFT,C.pharmaPCSAFT}
function x0_sat_pure_default(model,T) = @invoke C.x0_sat_pure(model::Any,T)

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

function can_superanc(model::SuperancPCSAFT)
    val = true
    val = val & (length(model.params.epsilon_assoc.values.values) != 0)
    if hasfield(typeof(model.params),:dipole)
        val = val & all(iszero,model.params.dipole.values)
    end
    return val && C.SUPERANC_ENABLED[]
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

x0_crit_pure_default(model) = @invoke C.x0_crit_pure(model::C.SAFTModel)

function C.x0_crit_pure(model::SuperancPCSAFT)
    can_superanc(model) || return x0_crit_pure_default(model)
    m,ϵ,σ = get_pcsaft_consts(model)
    if 1.0 <= m <= 64.0
        T̃c = pcsaft_tc(m,1.0)
        vc = pcsaft_vc(m,σ + Δσ(model,ϵ*T̃c))
        return T̃c,log10(vc)
    else
        return x0_crit_pure_default(model)
    end
end

function crit_pure(model::SuperancPCSAFT)
    can_superanc(model) || return crit_pure(model,C.x0_crit_pure(model))
    m,ϵ,σ = get_pcsaft_consts(model)
    if 1.0 <= m <= 64.0
        Tc = pcsaft_tc(m,ϵ)
        vc = pcsaft_vc(m,σ + Δσ(model,Tc))
        pc = pressure(model,vc,Tc)
        return Tc,pc,vc
    else
        return crit_pure(model,C.x0_crit_pure(model))
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

function can_superanc(model::SuperancCubic)
    return C.SUPERANC_ENABLED[]
end

function __init__()
    C.use_superancillaries!()
end

end #module
