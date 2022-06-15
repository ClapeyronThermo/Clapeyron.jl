
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
    !isone(length(model)) && throw(error("$model can only have one component."))
    T = T*T/T
    return saturation_pressure_impl(model,T,method)
end

include("ChemPotV.jl")
include("IsoFugacity.jl")
include("ChemPotDensity.jl")
include("SuperAnc.jl")
include("ClapeyronSat.jl")
include("AntoineSat.jl")

function saturation_temperature(model,p,method::SaturationMethod=ClapeyronSaturation{Nothing}(nothing))
    !isone(length(model)) && throw(error("$model can only have one component."))
    p = p*p/p
    return saturation_temperature_impl(model,p,method)
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


