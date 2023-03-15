struct PolExpLiquid <: LiquidVolumeModel
    Tc::Float64
    rhoc::Float64
    n::Vector{Float64}
    v::Vector{Float64}
end

function _rho_sat(model::PolExpLiquid,T)
    T_c = model.Tc
    ρ_c = model.rhoc
    T>T_c && return zero(T)/zero(T)
    Tr = T/T_c
    θ = 1.0-Tr
    ρ_r = evalexppoly(θ,model.n,model.v) + 1.0
    ρ_l = ρ_r*ρ_c
    return ρ_l
end

function volume_impl(model::PolExpLiquid,p,T,z::SingleComp,phase=:unknown,threaded=false,vol0 = 0.0)
    return 1/_rho_sat(model,T)
end

Base.length(::PolExpLiquid) = 1
is_splittable(::PolExpLiquid) = false

export PolExpLiquid