struct ClapeyronSaturation{T,C,M<:SaturationMethod} <: SaturationMethod
    T0::T
    crit::C
    satmethod::M
end

"""
    ClapeyronSaturation <: SaturationMethod
    ClapeyronSaturation(T0 = nothing, crit = nothing, satmethod = ChemPotVSaturation())

Saturation method for `saturation_temperature`. It solves iteratively `saturation_temperature(model,Ti,satmethod)` until convergence, by using the Clapeyron equation:
```
dp/dT = ΔS/ΔV
```
It descends from the critical point (or `T0`, if provided). Reliable, but slow.

It is recommended that `T0 > Tsat`, as the temperature decrease iteration series is more stable. Default method for `saturation_temperature` until Clapeyron 0.3.7
"""
ClapeyronSaturation

ClapeyronSaturation(;T0 = nothing, crit = nothing, satmethod = ChemPotVSaturation(;crit)) = ClapeyronSaturation(T0,crit,satmethod)
#if a model overloads x0_saturation_temperature to return a T0::Number, we can assume this number is near
#the actual saturation temperature, so we use the direct algorithm. otherwise, we use a safe approach, starting from the critical
#coordinate and descending.

#by default, starts right before the critical point, and descends via Clapeyron equation: (∂p/∂T)sat = ΔS/ΔV ≈ Δp/ΔT

function saturation_temperature_impl(model::EoSModel,p,method::ClapeyronSaturation{Nothing})
    crit = method.crit
    if isnothing(method.crit)
        crit = crit_pure(model)
    end
    SatMethod = parameterless_type(method.satmethod)
    satmethod = SatMethod(;crit)
    T0,_,_ = x0_saturation_temperature_crit(model,p,crit)
    method_init = ClapeyronSaturation(T0,crit,satmethod)
    return saturation_temperature_impl(model,p,method_init)
end

function saturation_temperature_impl(model::EoSModel,p,method::ClapeyronSaturation)
    T0 = method.T0/one(method.T0)
    TT = typeof(T0)
    nan = zero(T0)/zero(T0)
    cache = Ref{Tuple{TT,TT,TT,TT,Bool}}((nan,nan,nan,nan,false))
    f(T) = Obj_sat_pure_T(model,T,p,cache,method.satmethod)
    T = Solvers.fixpoint(f,T0)
    _,_,v_l,v_v,_ = cache[]
    return T,v_l,v_v
end

function Obj_sat_pure_T(model,T,p,cache,satmethod)
    Told,pold,_,_,_ = cache[]
    pii,vli,vvi = saturation_pressure(model,T,satmethod)
    Ti,sat = refine_x0_saturation_temperature(model,p,T,satmethod,false)
    pii,vli,vvi = sat
    Δp = (p-pii)
    abs(Δp) < 4eps(p) && return T
    #if abs(Δp/p) < 0.01
    #    use_v = true
    #end
    if Told < T
        if isnan(pii) && !isnan(pold)
            return (T+Told)/2
        end
    end
    cache[] = (T,pii,vli,vvi,false)
    return Ti
end

export ClapeyronSaturation