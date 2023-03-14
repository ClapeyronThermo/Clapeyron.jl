struct PolExpVapour <: LiquidVolumeModel
    Tc::Float64
    rhoc::Float64
    n::Vector{Float64}
    v::Vector{Float64}
end

function _rho_sat(model::PolExpVapour,T)
    T_c = model.Tc
    ρ_c = model.rhoc
    T>T_c && return zero(T)/zero(T)
    Tr = T/T_c
    θ = 1.0-Tr
    log_ρ_v_ρ_c =evalexppoly(θ,model.n,model.v)
    ρ_v = exp(log_ρ_v_ρ_c)*ρ_c
    return ρ_v
end

function volume_impl(model::PolExpVapour,p,T,z::SingleComp,phase=:unknown,threaded=false,vol0 = 0.0)
    return 1/_rho_sat(model,T)
end

Base.length(::PolExpVapour) = 1
is_splittable(::PolExpVapour) = false