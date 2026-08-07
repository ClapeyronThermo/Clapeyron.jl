struct PolExpLiquid <: LiquidVolumeModel
    data::GenericAncEvaluator
end

PolExpLiquid(anc::Solvers.ChebyshevRange) = PolExpLiquid(GenericAncEvaluator(anc))

function _rho_sat(model::PolExpLiquid,T)
    return _eval_generic_anc(model.data,T)
end

function volume_impl(model::PolExpLiquid,p,T,z,phase,threaded,vol0)
    @assert length(z) == 1
    return z[1]/_rho_sat(model,T)
end

Base.length(::PolExpLiquid) = 1
is_splittable(::PolExpLiquid) = false

export PolExpLiquid