#Tillner-Roth and Friend mixing rule
struct TillnerRothFriendMixing <: MixingRule end
TillnerRothFriendMixing(::Any, ::Any) = TillnerRothFriendMixing()

# IAPWS-95 water Tc = 647.096 K; ammonia Tc ≈ 405.4 K → threshold 500 K is safe.
@inline _is_water(model::MultiFluid, i::Int) = model.pures[i].properties.Tc > 500.0

@inline function _water_ammonia_indices(model::MultiFluid)
    iw = _is_water(model, 1) ? 1 : 2
    ia = 3 - iw
    return iw, ia
end

function T_scale(model::MultiFluid, z, ::TillnerRothFriendMixing, ∑z)
    Tc = model.params.Tr.values
    length(model) == 1 && return only(Tc)*one(∑z)
    kt = 0.9648407
    α = 1.125455
    iw, ia = _water_ammonia_indices(model)

    Tcw, Tca = Tc[iw], Tc[ia]
    xw, xa = z[iw]/∑z, z[ia]/∑z
    T12_2 = kt*(Tcw + Tca) #2*T12
    return xw*xw*Tcw + xa*xa*Tca + xa * (1 - xa^α) * T12_2
end

function v_scale(model::MultiFluid, z, ::TillnerRothFriendMixing, ∑z)
    Vc = model.params.Vr.values
    length(model) == 1 && return only(Vc)*one(∑z)
    kv = 1.2395117
    β = 0.8978069
    iw, ia = _water_ammonia_indices(model)
    Vcw, Vca = Vc[iw], Vc[ia]
    xw, xa = z[iw]/∑z, z[ia]/∑z
    V12_2 = kv*(Vcw + Vca) #2*V12
    return xw*xw*Vcw + xa*xa*Vca + xa * (1 - xa^β) * V12_2
end

is_splittable(::TillnerRothFriendMixing) = false

export TillnerRothFriendMixing
