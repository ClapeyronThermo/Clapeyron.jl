
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
    T = T*(T/T)*oneunit(eltype(model))
    satmodel = saturation_model(model)
    satmodel !== model && saturation_pressure(satmodel,T,method)
    if has_a_res(model)
        res = saturation_pressure_impl(primalval(model),primalval(T),method)
        saturation_pressure_ad(model,T,res)
    else
        return saturation_pressure_impl(model,T,method)
    end
end

function saturation_pressure(model::EoSModel,T;kwargs...)
    satmodel = saturation_model(model)
    if satmodel !== model
        return saturation_pressure(satmodel,T;kwargs...)
    end
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
    satmodel = saturation_model(model)
    if satmodel !== model
        return saturation_pressure(satmodel,T,V0)
    end
    single_component_check(saturation_pressure,model)
    vl = first(V0)
    vv = last(V0)
    kwargs = (;vl,vv)
    method = init_preferred_method(saturation_pressure,model,kwargs)
    return saturation_pressure(model,T,method)
end

function saturation_pressure_ad(model,T,result)
    if has_dual(model) || has_dual(T)
        p_primal,vl_primal,vv_primal = result

        #=
        update step from https://github.com/lucpaoli/SAFT_ML/blob/f22648055bdf4cd244cf427a596fc7b1c03e6383/saftvrmienn.jl#L138-L157
        =#
        Δg = eos(model, vv_primal, T) - eos(model,vl_primal, T) + p_primal*(vv_primal - vl_primal)
        Δv = vv_primal - vl_primal
        p = p_primal - Δg/Δv
        #=
        for volume, we use a volume update
        =#
        vl = volume_ad(model,vl_primal,T,SA[1.0],p)
        vv = volume_ad(model,vv_primal,T,SA[1.0],p)
        return p,vl,vv
    else
        return result
    end
end

function derivx(f,i)
    h = sqrt(eps(i))
    return (1*f(i-2h)-8*f(i-1h)+0*f(i)+8*f(i+1h)-1*f(i+2h))/(12*h)
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
    satmodel = saturation_model(model)
    if satmodel !== model
        return saturation_temperature(satmodel,p;kwargs...)
    end
    method = init_preferred_method(saturation_temperature,model,kwargs)
    return saturation_temperature(model,p,method)
end

function saturation_temperature(model,p,method::SaturationMethod)
    satmodel = saturation_model(model)
    if satmodel !== model
        return saturation_temperature(satmodel,p;method)
    end
    single_component_check(crit_pure,model)
    p = p*p/p

    if has_a_res(model)
        res = saturation_temperature_impl(model,primalval(p),method)
        return saturation_temperature_ad(model,p,res)
    else
        return saturation_temperature_impl(model,p,method)
    end
end

#if a number is provided as initial point, it will instead proceed to solve directly
function saturation_temperature(model::EoSModel, p, T0::Number)
    satmodel = saturation_model(model)
    if satmodel !== model
        return saturation_temperature(satmodel,p,T0)
    end
    kwargs = (;T0)
    method = init_preferred_method(saturation_temperature,model,kwargs)
    saturation_temperature(model,p,method)
end

function saturation_temperature_ad(model,p,result)
    if has_dual(model) || has_dual(p)
        T_primal,vl_primal,vv_primal = result
        vl = volume_ad(model,vl_primal,T_primal,SA[1.0],p)

        #manual volume_ad for vapour volume, we reuse p_primal for the calculation of the T step.
        p_primal,∂p∂V = p∂p∂V(model,vv_primal,T_primal,SA[1.0])
        vv = vv_primal - (p_primal - p)/∂p∂V

        #for T, we use a dlnpdTinv step, a dpdT step is fine too
        dpdT = dpdT_saturation(model,vv_primal,vl_primal,T_primal)
        dTinvdlnp = -p_primal/(dpdT*T_primal*T_primal)
        Δlnp = log(p/p_primal)
        Tinv0 = 1/T_primal
        Tinv = Tinv0 + dTinvdlnp*Δlnp
        dT = T_primal - 1/Tinv
        T = 1/Tinv
        T = T_primal - (p_primal - p)/dpdT
        return T,vl,vv
    else
        return result
    end
end

include("ChemPotV.jl")
include("IsoFugacity.jl")
include("ChemPotDensity.jl")
include("SuperAnc.jl")
include("CritExtrapolation.jl")
include("ClapeyronSat.jl")
include("AntoineSat.jl")


"""
    enthalpy_vap(model::EoSModel, T,method = ChemPotVSaturation(x0_sat_pure(model,T)))

Calculates `ΔH`, the difference between saturated vapour and liquid enthalpies at temperature `T`, in J
"""
function enthalpy_vap(model::EoSModel, T,satmethod = ChemPotVSaturation())
    single_component_check(enthalpy_vap,model)
    (P_sat,V_l,V_v) = saturation_pressure(model,T,satmethod)
    H_v = VT_enthalpy_res(model,V_v,T)
    H_l = VT_enthalpy_res(model,V_l,T)
    #H_v(res) - H_l(res) = H_l - H_v
    H_vap = H_v - H_l
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

function fast_pcsaft1(eps::T,sigma::T,segment::T,mw::T = oneunit(T)) where T
    name1 = ""
    empty_str = String[]
    namevec = [name1]
    
    missingval_vec = [false]
    missingval_mat = [false;;]
    empty_nested_int = [Int64[]]
    empty_nested_string = [empty_str]
    assoc_vec = Compressed4DMatrix{T}()
    ϵ = PairParameter(name1,namevec,[eps;;],missingval_mat,empty_str,empty_str)
    σ = PairParameter(name1,namevec,[sigma;;],missingval_mat,empty_str,empty_str)
    m = SingleParameter(name1,namevec,[segment],missingval_vec,empty_str,empty_str)
    Mw = SingleParameter(name1,namevec,[mw],missingval_vec,empty_str,empty_str)
    bondvol = AssocParam(name1,namevec,assoc_vec,empty_nested_string,empty_str,empty_str)
    ϵ_ab = AssocParam(name1,namevec,assoc_vec,empty_nested_string,empty_str,empty_str)

    param = PCSAFTParam{T}(Mw,m,σ,ϵ,ϵ_ab,bondvol)
    
    sites = SiteParam(namevec,[empty_str],PackedVectorsOfVectors.pack(empty_nested_int),empty_nested_int,empty_str,empty_nested_int,empty_str,empty_str,nothing)
    return PCSAFT{BasicIdeal,T}(namevec,sites,param,BasicIdeal(),AssocOptions(),empty_str)
end

function test_f(x)
    eps,sigma,segment = x[1],x[2],x[3]
    model = fast_pcsaft1(100*eps,sigma*1e-10,segment)
    p,vl,vv = saturation_pressure(model,140.15)
    return p
end
