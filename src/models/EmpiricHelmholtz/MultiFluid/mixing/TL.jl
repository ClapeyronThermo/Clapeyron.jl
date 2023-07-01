#Tillner-Roth and Friend mixing rule
struct TillnerRothFriendMixing <: MixingRule
    components::Vector{String}
    nh3_idx::SpecialComp
end

function v_scale(model::MultiFluid,z,mixing::TillnerRothFriendMixing,∑z)
    Vc = model.params.Vc.values
    length(model) == 1 && return only(Vc)
    ∑z⁻¹ = 1/∑z
    if mixing.nh3_idx[] == 1
        xa,xw = z[1]*∑z⁻¹,z[2]*∑z⁻¹
        Vca,Vcw = Vc[1],Vc[2]
    else
        xa,xw = z[2]*∑z⁻¹,z[1]*∑z⁻¹
        Vca,Vcw = Vc[2],Vc[1]
    end
    kv = 1.2395117
    β = 0.8978069
    return xw*xw*Vcw + xa*xa*Vca + kv*(1-xw)*(1-xa^β)*(Vca + Vcw)
end

function T_scale(model::MultiFluid,z,mixing::TillnerRothFriendMixing,∑z)
    Tc = model.params.Tc.values
    length(model) == 1 && return only(Tc)
    ∑z⁻¹ = 1/∑z
    if mixing.nh3_idx[] == 1
        xa,xw = z[1]*∑z⁻¹,z[2]*∑z⁻¹
        Tca,Tcw = Tc[1],Tc[2]
    else
        xa,xw = z[2]*∑z⁻¹,z[1]*∑z⁻¹
        Tca,Tcw = Tc[2],Tc[1]
    end
    kt = 0.9648407
    α = 1.125455
    return xw*xw*Tcw + xa*xa*Tca + kt*(1-xw)*(1-xa^α)*(Tca + Tcw)
end

export TillnerRothFriendMixing