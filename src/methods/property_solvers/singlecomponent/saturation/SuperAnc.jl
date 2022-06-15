
struct SuperAncSaturation <: SaturationMethod 
    p_tol::Float64
end

function SuperAncSaturation(;p_tol = 1e-16)
    return SuperAncSaturation(p_tol)
end

function saturation_pressure_impl(model::ABCubicModel,T,method::SuperAncSaturation)
    Tc = model.params.Tc.values[1]
    a,b,c = cubic_ab(model,1e-3,T)
    T̃ = T*R̄*b/a

    Vv = chebyshev_vapour_volume(model,T̃,b)
    Vl = chebyshev_liquid_volume(model,T̃,b)
    p_sat = chebyshev_pressure(model,T̃,a,b)
    return (p_sat,Vl,Vv)
end

function chebyshev_vapour_volume(model::ABCubicModel,T̃,b)
    Cₙ = chebyshev_coef_v(model)
    Tmin = chebyshev_Tmin_v(model)
    Tmax = chebyshev_Tmax_v(model)
    for i ∈ 1:length(Cₙ)
        if T̃<=Tmax[i] && T̃>=Tmin[i]
            c = ChebyshevT(Cₙ[i])
            T̄ = (2*T̃ - (Tmax[i] + Tmin[i])) / (Tmax[i] - Tmin[i])
            ρ̃ᵥ = c(T̄)            
            Vᵥ = b/ρ̃ᵥ
            return Vᵥ
        end
    end
end

function chebyshev_liquid_volume(model::ABCubicModel,T̃,b)
    Cₙ = chebyshev_coef_l(model)
    Tmin = chebyshev_Tmin_l(model)
    Tmax = chebyshev_Tmax_l(model)
    for i ∈ 1:length(Cₙ)
        if T̃<=Tmax[i] && T̃>=Tmin[i]
            c = ChebyshevT(Cₙ[i])
            T̄ = (2*T̃ - (Tmax[i] + Tmin[i])) / (Tmax[i] - Tmin[i])
            ρ̃ₗ = c(T̄)
            Vₗ = b/ρ̃ₗ
            return Vₗ
        end
    end
end

function chebyshev_pressure(model::ABCubicModel,T̃,a,b)
    Cₙ = chebyshev_coef_p(model)
    Tmin = chebyshev_Tmin_p(model)
    Tmax = chebyshev_Tmax_p(model)
    for i ∈ 1:length(Cₙ)
        if T̃<=Tmax[i] && T̃>=Tmin[i]
            c = ChebyshevT(Cₙ[i])
            T̄ = (2*T̃ - (Tmax[i] + Tmin[i])) / (Tmax[i] - Tmin[i])
            p̃ = c(T̄)
            p = p̃*a/b^2
            return p
        end
    end
end

export SuperAncSaturation
