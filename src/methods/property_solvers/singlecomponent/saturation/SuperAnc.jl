
struct SuperAncSaturation <: SaturationMethod
    p_tol::Float64
end

const SUPERANC_ENABLED = Ref(false)

"""
    use_superancillaries!(val::Bool = true)

Enable the use of cubic and PC-SAFT superancillaries as initial points for `saturation_pressure`. For `PCSAFT` it also enables the use of superancillaries for critical point calculations.
This function requires `EoSSuperancillaries.jl` to be loaded in the current session (`using EoSSuperancillaries`).
"""
function use_superancillaries!(val::Bool = true)
    SUPERANC_ENABLED[] = val
    return val
end

"""
    SuperAncSaturation <: SaturationMethod
    SuperAncSaturation()

Saturation method for `saturation_pressure`. It uses Chebyshev expansions constructed with extended precision arithmetic to obtain the saturation curves of supported EoS within `Float64` precision. Models supported are:
 - [vdW](@ref) models and variants
 - [RK](@ref) models and variants
 - [PR](@ref) models and variants
 - [PCSAFT](@ref) models and some variants (via `EoSSuperancillaries.jl` package extension)

## References
1. Bell, I. H., & Deiters, U. K. (2021). Superancillary equations for cubic equations of state. Industrial & Engineering Chemistry Research, 60(27), 9983–9991. doi:10.1021/acs.iecr.1c00847
"""
SuperAncSaturation

function SuperAncSaturation(;p_tol = 1e-16,crit = nothing)
    return SuperAncSaturation(p_tol)
end

#=
TODO: move everything to the extension on the next breaking change
=#

function saturation_pressure_impl(model::ABCubicModel,T,method::SuperAncSaturation)
    Tc = model.params.Tc.values[1]
    if Tc < T
        nan = zero(T)/zero(T)
        return (nan,nan,nan)
    end
    a,b,c = cubic_ab(model,1e-3,T)
    T̃ = T*Rgas(model)*b/a
    Vv = chebyshev_vapour_volume(model,T̃,b) - c
    Vl = chebyshev_liquid_volume(model,T̃,b) - c
    p_sat = chebyshev_pressure(model,T̃,a,b)
    return (p_sat,Vl,Vv)
end

function chebyshev_vapour_volume(model::ABCubicModel,T̃,b)
    cheb_vvsat = _cheb_vvsat(model)
    return b/Solvers.cheb_eval(cheb_vvsat,T̃)
end

function chebyshev_liquid_volume(model::ABCubicModel,T̃,b)
    cheb_vlsat = _cheb_vlsat(model)
    return b/Solvers.cheb_eval(cheb_vlsat,T̃)
end

function chebyshev_pressure(model::ABCubicModel,T̃,a,b)
    cheb_psat = _cheb_psat(model)
    p̃ = Solvers.cheb_eval(cheb_psat,T̃)
    return p̃*a/(b*b)
end

function saturation_temperature_impl(model::ABCubicModel,p,method::SuperAncSaturation)
    Tc = model.params.Tc.values[1]
    Pc = model.params.Pc.values[1]
    nan = zero(p)/zero(p)
    if Pc < p
        return (nan,nan,nan)
    end
    T = chebyshev_temperature(model,p,method)
    a,b,c = cubic_ab(model,1e-3,T)
    T̃ = T*Rgas(model)*b/a
    Vv = chebyshev_vapour_volume(model,T̃,b) - c
    Vl = chebyshev_liquid_volume(model,T̃,b) - c
    return (T,Vl,Vv)
    #p_sat = chebyshev_pressure(model,T̃,a,b)
end

function chebyshev_temperature(model::ABCubicModel,p,method::SuperAncSaturation)
    
    function f0(τ)
        T = 1/τ
        a,b,c = cubic_ab(model,1e-3,T)
        T̃ = T*Rgas(model)*b/a
        return chebyshev_pressure(model,T̃,a,b) - p
    end

    Tc,Pc,Vc = crit_pure(model)
    a07,b07,_ = cubic_ab(model,1e-3,0.7Tc)
    T̃07 = 0.7*Tc*Rgas(model)*b07/a07
    p07 = chebyshev_pressure(model,T̃07,a07,b07)
    ω = -log10(p07/Pc) - 1.0
    T0 = Tc/(1 - log(p/Pc)/(5.3726985503194395*(1+ω)))
    prob = Roots.ZeroProblem(f0,1/T0)
    T = 1/Roots.solve(prob)
end

export SuperAncSaturation, use_superancillaries!
