
"""
    SaturationMethod

Abstract type for `saturation_temperature` and `saturation_pressure` routines.

"""

abstract type SaturationMethod end


"""
    saturation_pressure(model::EoSModel, T)
    saturation_pressure(model::EoSModel,T,method::SaturationMethod)
    saturation_pressure(model,T,x0::Union{Tuple,Vector})

Performs a single component saturation equilibrium calculation, at the specified temperature `T`, of one mol of pure sustance specified by `model`

Returns `(p₀, Vₗ, Vᵥ)` where `p₀` is the saturation pressure (in Pa), `Vₗ` is the liquid saturation volume (in m³) and `Vᵥ` is the vapour saturation volume (in m³).

If the calculation fails, returns  `(NaN, NaN, NaN)`

By default, it uses [`ChemPotVSaturation`](@ref)
## Examples:
julia-repl
```
julia> pr = PR(["water"])
PR{BasicIdeal, PRAlpha, NoTranslation, vdW1fRule} with 1 component:
 "water"
Contains parameters: a, b, Tc, Pc, Mw

julia> saturation_pressure(pr,373.15) #default, uses Clapeyron.ChemPotVSaturation
(96099.38979351855, 2.2674781912892906e-5, 0.03201681565699426)

julia> saturation_pressure(pr,373.15,IsoFugacitySaturation()) #iso fugacity
(96099.38979351871, 2.2674781912892933e-5, 0.03201681565699359)

julia> saturation_pressure(pr,373.15,IsoFugacitySaturation(p0 = 1.0e5)) #iso fugacity, with starting point
(96099.38979351871, 2.2674781912892933e-5, 0.03201681565699547)
```

"""

function saturation_pressure(model,T,method::SaturationMethod)
    !isone(length(model)) && throw(error("$model have more than one component."))
    T = T*T/T
    return saturation_pressure_impl(model,T,method)
end

include("ChemPotV.jl")
include("IsoFugacity.jl")
include("ChemPotDensity.jl")

#if a model overloads x0_saturation_temperature to return a T0::Number, we can assume this number is near
#the actual saturation temperature, so we use the direct algorithm. otherwise, we use a safe approach, starting from the critical
#coordinate and descending.
x0_saturation_temperature(model,p) = nothing

#by default, starts right before the critical point, and descends via Clapeyron equation: (∂p/∂T)sat = ΔS/ΔV ≈ Δp/ΔT
saturation_temperature(model::EoSModel, P) = saturation_temperature(model, P,x0_saturation_temperature(model,P))

function saturation_temperature(model::EoSModel,p,T0::Nothing)
    Tc,pc,vc = crit_pure(model)
    TT = typeof(vc*pc/Tc)
    nan = zero(TT)/zero(TT)
    p > 0.99999pc && (return (nan,nan,nan)) 
    T0 = 0.99*Tc
    isnan(T0) && (return (nan,nan,nan))
    cache = Ref{Tuple{TT,TT,TT,TT}}((nan,nan,nan,nan))
    f(T) = Obj_sat_pure_T(model,T,p,cache)
    T = Solvers.fixpoint(f,T0)
    _,_,v_l,v_v = cache[]
    return T,v_l,v_v
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

function obj_tsat(model::EoSModel, T, P,cache)
    vol_liq,vol_vap = cache[]
    RT = R̄*T

    vol_liq = volume(model, P, T, phase=:liquid, vol0=vol_liq)
    vol_vap = volume(model, P, T, phase=:vapor, vol0=vol_vap)

    μ_liq = VT_chemical_potential_res(model, vol_liq, T)[1]
    μ_vap = VT_chemical_potential_res(model, vol_vap, T)[1]

    Z_liq = P*vol_liq/RT
    Z_vap = P*vol_vap/RT

    lnϕ_liq = μ_liq/RT - log(Z_liq)
    lnϕ_vap = μ_vap/RT - log(Z_vap)
    FO = lnϕ_vap - lnϕ_liq
    cache[] = (vol_liq,vol_vap)
    return FO
end

#if a number is provided as initial point, it will instead proceed to solve directly
function saturation_temperature(model::EoSModel, P, T0::Number)
    vol_liq = volume(model, P, T0, phase=:liquid)
    vol_vap = volume(model, P, T0, phase=:vapor)
    cache = Ref((vol_liq,vol_vap))
    ftsat(T) = obj_tsat(model, T, P,cache)
    fT = Roots.ZeroProblem(ftsat,T0)
    T = Roots.solve(fT,Roots.Order0())
    vol_liq,vol_vap = cache[]
    return T, vol_liq, vol_vap
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

#tsat, psat interface
include("tsat_psat.jl")


