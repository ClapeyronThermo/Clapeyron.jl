struct PolExpVapour <: LiquidVolumeModel
    data::GenericAncEvaluator
end

function _rho_sat(model::PolExpVapour,T)
    return _eval_generic_anc(model.data,T)
end

function volume_impl(model::PolExpVapour,p,T,z,phase,threaded,vol0)
    @assert length(z) == 1
    return z[1]/_rho_sat(model,T)
end

Base.length(::PolExpVapour) = 1
is_splittable(::PolExpVapour) = false

export PolExpVapour