
"""
    SaturationMethod

Abstract type for `saturation_temperature` and `saturation_pressure` routines.

"""

abstract type SaturationMethod end


"""
    saturation_pressure(model::EoSModel, T, V0 = x0_sat_pure(model,T))

Performs a single component saturation equilibrium calculation, at the specified temperature `T`, of one mol of pure sustance specified by `model`

Returns `(p₀, Vₗ, Vᵥ)` where `p₀` is the saturation pressure (in Pa), `Vₗ` is the liquid saturation volume (in m³) and `Vᵥ` is the vapour saturation volume (in m³).

If the calculation fails, returns  `(NaN, NaN, NaN)`

`V0` is `[log10(Vₗ₀),log10(Vᵥ₀)]` , where `Vₗ₀`  and `Vᵥ₀` are initial guesses for the liquid and vapour volumes.
"""

function saturation_pressure(model,T,method::SaturationMethod)
    !isone(length(model)) && throw(error("$model have more than one component."))
    T = T*T/T
    return saturation_pressure_impl(model,T,method)
end

include("ChemPotV.jl")

#by default, starts right before the critical point, and descends via Clapeyron equation: (∂p/∂T)sat = ΔS/ΔV ≈ Δp/ΔT

function saturation_temperature(model::EoSModel,p,T0 = nothing)
    if T0 === nothing
        T0 = x0_saturation_temperature(model,p)
    end
    TT = typeof(T0/T0)
    nan = zero(TT)/zero(TT)
    if isnan(T0)
        return (nan,nan,nan)
    end 
    cache = Ref{Tuple{TT,TT,TT,TT}}((nan,nan,nan,nan))
    f(T) = Obj_sat_pure_T(model,T,p,cache)
    T = Solvers.fixpoint(f,T0)
    _,_,v_l,v_v = cache[]
    return T,v_l,v_v
end

function x0_saturation_temperature(model,p)
    Tc,pc,vc = crit_pure(model)    
    p > 0.99999pc && (return zero(p)/zero(p))
    T99 = 0.99*Tc
    return T99
end

function Obj_sat_pure_T(model,T,p,cache)
    Told,pold,vlold,vvold = cache[]
    pii,vli,vvi = saturation_pressure(model,T)
    Δp = (p-pii)
    abs(Δp) < 4eps(p) && return T
    if Told < T
        if isnan(pii) && !isnan(pold)
            return (T+Told)/2
        end
    end
    cache[] = (T,pii,vli,vvi) 
    S_v = VT_entropy(model,vvi,T)
    S_l = VT_entropy(model,vli,T)
    ΔS = S_v - S_l
    ΔV = vvi - vli
    dpdt = ΔS/ΔV #≈ (p - pii)/(T-Tnew)
    Ti = T + Δp/dpdt
    return Ti
end

"""
    enthalpy_vap(model::EoSModel, T,method = ChemPotVSaturation(x0_sat_pure(model,T)))

Calculates `ΔH`, the difference between saturated vapour and liquid enthalpies at temperature `T`, in J   
"""
function enthalpy_vap(model::EoSModel, T,method = ChemPotVSaturation(x0_sat_pure(model,T)))
    (P_sat,V_l,V_v) = saturation_pressure(model,T,method)
    H_v = VT_enthalpy(model,V_v,T)
    H_l = VT_enthalpy(model,V_l,T)
    H_vap=H_v -H_l
    return H_vap
end

#TODO: support method as optional parameter
"""
    acentric_factor(model::EoSModel)

calculates the acentric factor using its definition:

    ω = -log10(psatᵣ) -1, at Tᵣ = 0.7
To do so, it calculates the critical temperature (using `crit_pure`) and performs a saturation calculation (with `sat_pure`)

"""
function acentric_factor(model::EoSModel)
    T_c,p_c,_ = crit_pure(model)
    p = first(saturation_pressure(model,0.7*T_c))
    p_r = p/p_c
    return -log10(p_r) - 1.0
end



