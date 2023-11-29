
"""
    SaturationMethod <: ThermodynamicMethod 

Abstract type for `saturation_temperature` and `saturation_pressure` routines.
Should at least support passing the `crit` keyword, containing the critical point, if available.
"""
abstract type SaturationMethod <: ThermodynamicMethod end

"""
    saturation_pressure(model::EoSModel, T)
    saturation_pressure(model::EoSModel,T,method::SaturationMethod)
    saturation_pressure(model,T,x0::Union{Tuple,Vector})

Performs a single component saturation equilibrium calculation, at the specified temperature `T`, of one mol of pure sustance specified by `model`
Returns `(p₀, Vₗ, Vᵥ)` where `p₀` is the saturation pressure (in Pa), `Vₗ` is the liquid saturation volume (in m³) and `Vᵥ` is the vapour saturation volume (in m³).

If the calculation fails, returns  `(NaN, NaN, NaN)`

By default, it uses [`ChemPotVSaturation`](@ref)

## Examples:

```julia-repl
julia> pr = PR(["water"])
PR{BasicIdeal, PRAlpha, NoTranslation, vdW1fRule} with 1 component:
  "water"
Contains parameters: a, b, Tc, Pc, Mw

julia> p,vl,vv = saturation_pressure(pr,373.15) #default, uses Clapeyron.ChemPotVSaturation
(96099.38979351855, 2.2674781912892906e-5, 0.03201681565699426)

julia> p,vl,vv = saturation_pressure(pr,373.15,IsoFugacitySaturation()) #iso fugacity
(96099.38979351871, 2.2674781912892933e-5, 0.03201681565699359)

julia> p,vl,vv = saturation_pressure(pr,373.15,IsoFugacitySaturation(p0 = 1.0e5)) #iso fugacity, with starting point
(96099.38979351871, 2.2674781912892933e-5, 0.03201681565699547)
```
"""
function saturation_pressure(model::EoSModel,T,method::SaturationMethod)
    single_component_check(saturation_pressure,model)
    T = T*(T/T)
    return saturation_pressure_impl(model,T,method)
end

function saturation_pressure(model::EoSModel,T;kwargs...)
    if keys(kwargs) == (:v0,)
        nt_kwargs = NamedTuple(kwargs)
        v0 = nt_kwargs.v0
        vl = first(v0)
        vv = last(v0)
        _kwargs = (;vl,vv)
        method = init_preferred_method(saturation_pressure,model,_kwargs)
    else
        method = init_preferred_method(saturation_pressure,model,kwargs)
    end
    return saturation_pressure(model,T,method)
end

function saturation_pressure(model::EoSModel,T,V0::Union{Tuple,Vector})
    single_component_check(saturation_pressure,model)
    vl = first(V0)
    vv = last(V0)
    kwargs = (;vl,vv)
    method = init_preferred_method(saturation_pressure,model,kwargs)
    return saturation_pressure(model,T,method)
end


"""
    check_valid_sat_pure(model,P_sat,Vl,Vv,T,ε0 = 5e7)

Checks that a saturation method converged correctly. it checks:
- That both volumes are mechanically stable
- That both volumes are different, with a difference of at least `ε0` epsilons
"""
function check_valid_sat_pure(model,P_sat,V_l,V_v,T,ε0 = 5e7)
    ε = abs(V_l-V_v)/(eps(typeof(V_l-V_v)))
    ε <= ε0 && return false
    pl,dpdvl = p∂p∂V(model,V_l,T,SA[1.0])
    pv,dpdvv = p∂p∂V(model,V_v,T,SA[1.0])
    ε1 = ε0*eps(P_sat)
    return  (dpdvl <= 0)        && #mechanical stability of the liquid phase
            (dpdvv <= 0)        && #mechanical stability of the vapour phase
            (1 - pl/P_sat < ε1)  && #calculated pressure for liquid phase is aproximately the same
            (1 - pv/P_sat < ε1)     #calculated pressure for vapour phase is aproximately the same
    #if ΔV > ε then Vl and Vv are different values
end

"""
    saturation_temperature(model::EoSModel, p, kwargs...)
    saturation_temperature(model::EoSModel, p, method::SaturationMethod)
    saturation_temperature(model, p, T0::Number)

Performs a single component saturation temperature equilibrium calculation, at the specified pressure `T`, of one mol of pure sustance specified by `model`

Returns `(T₀, Vₗ, Vᵥ)` where `p₀` is the saturation Temperature (in K), `Vₗ` is the liquid saturation volume (in m³) and `Vᵥ` is the vapour saturation volume (in m³).

If the calculation fails, returns  `(NaN, NaN, NaN)`

By default, it uses [`AntoineSaturation`](@ref)

## Examples:
julia-repl
```
julia> pr = PR(["water"])
PR{BasicIdeal, PRAlpha, NoTranslation, vdW1fRule} with 1 component:
  "water"
Contains parameters: a, b, Tc, Pc, Mw

julia> Ts,vl,vv = saturation_temperature(pr,1e5) # AntoineSaturation by default
(374.24014010712983, 2.269760164801948e-5, 0.030849387955737825)

julia> saturation_pressure(pr,Ts)
(100000.00004314569, 2.269760164804427e-5, 0.03084938795785433)
```
"""
function saturation_temperature(model,p;kwargs...)
    method = init_preferred_method(saturation_temperature,model,kwargs)
    return saturation_temperature(model,p,method)
end

function saturation_temperature(model,p,method::SaturationMethod)
    single_component_check(crit_pure,model)
    p = p*p/p
    return saturation_temperature_impl(model,p,method)
end

#if a number is provided as initial point, it will instead proceed to solve directly
function saturation_temperature(model::EoSModel, p, T0::Number)
    kwargs = (;T0)
    method = init_preferred_method(saturation_temperature,model,kwargs)
    saturation_temperature(model,p,method)
end

include("ChemPotV.jl")
include("IsoFugacity.jl")
include("ChemPotDensity.jl")
include("SuperAnc.jl")
include("ClapeyronSat.jl")
include("AntoineSat.jl")


"""
    enthalpy_vap(model::EoSModel, T,method = ChemPotVSaturation(x0_sat_pure(model,T)))

Calculates `ΔH`, the difference between saturated vapour and liquid enthalpies at temperature `T`, in J   
"""
function enthalpy_vap(model::EoSModel, T,satmethod = ChemPotVSaturation())
    single_component_check(enthalpy_vap,model)
    (P_sat,V_l,V_v) = saturation_pressure(model,T,satmethod)
    H_v = VT_enthalpy(model,V_v,T)
    H_l = VT_enthalpy(model,V_l,T)
    H_vap=H_v -H_l
    return H_vap
end

"""
    acentric_factor(model::EoSModel;crit = crit_pure(model), satmethod = ChemPotVSaturation())

calculates the acentric factor using its definition:
```
ω = -log10(psatᵣ) -1, at Tᵣ = 0.7
```

To do so, it calculates the critical temperature (using `crit_pure`) and performs a saturation calculation (with `saturation_pressure(model,0.7Tc,satmethod)`)
"""
function acentric_factor(model::EoSModel;crit = crit_pure(model),satmethod = ChemPotVSaturation())
    return acentric_factor(model,crit,satmethod)
end

function acentric_factor(model::EoSModel,crit,satmethod)
    single_component_check(acentric_factor,model)
    T_c,p_c,_ = crit
    p = first(saturation_pressure(model,0.7*T_c,satmethod))
    p_r = p/p_c
    return -log10(p_r) - 1.0
end

function saturation_liquid_density(model::EoSModel,T,satmethod = ChemPotVSaturation())
    single_component_check(saturation_liquid_density,model)
    return saturation_pressure(model,T,satmethod)[2]
end

#default initializers for saturation pressure and saturation temperature

function init_preferred_method(method::typeof(saturation_pressure),model::EoSModel,kwargs)
    ChemPotVSaturation(;kwargs...)
end

function init_preferred_method(method::typeof(saturation_temperature),model::EoSModel,kwargs)
    return AntoineSaturation(;kwargs...)
end
