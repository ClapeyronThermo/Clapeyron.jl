#Tillner-Roth and Friend mixing rule
struct TillnerRothFriendMixing <: MixingRule end
TillnerRothFriendMixing(::Any, ::Any) = TillnerRothFriendMixing()

# IAPWS-95 water Tc = 647.096 K; ammonia Tc ≈ 405.4 K → threshold 500 K is safe.
@inline _is_water(model::MultiFluid, i::Int) =
    model.pures[i].properties.Tc > 500.0

@inline function _water_ammonia_indices(model::MultiFluid)
    iw = _is_water(model, 1) ? 1 : 2
    ia = 3 - iw
    return iw, ia
end

function T_scale(model::MultiFluid, z, ::TillnerRothFriendMixing, sumz)
    k_T   = 0.9648407
    alpha = 1.125455
    iw, ia = _water_ammonia_indices(model)
    Tc_w = model.pures[iw].properties.Tc
    Tc_a = model.pures[ia].properties.Tc
    xw, xa = z[iw], z[ia]
    T12 = (1.0 - k_T * xa^alpha) * 0.5 * (Tc_w + Tc_a)
    return (xw^2 * Tc_w + xa^2 * Tc_a + 2.0 * xw * xa * T12) / sumz^2
end

function v_scale(model::MultiFluid, z, ::TillnerRothFriendMixing, sumz)
    k_V  = 1.2395117
    beta = 0.8978069
    iw, ia = _water_ammonia_indices(model)
    Vc_w = 1/model.pures[iw].properties.rhoc
    Vc_a = 1/model.pures[ia].properties.rhoc
    xw, xa = z[iw], z[ia]
    V12 = (1.0 - k_V * xa^beta) * 0.5 * (Vc_w + Vc_a)
    return (xw^2 * Vc_w + xa^2 * Vc_a + 2.0 * xw * xa * V12) / sumz^2
end

export TillnerRothFriendMixing