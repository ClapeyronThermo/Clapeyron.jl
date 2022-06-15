struct ClapeyronSaturation{T} <: SaturationMethod
    Temp::Union{Nothing,T}
end

#if a model overloads x0_saturation_temperature to return a T0::Number, we can assume this number is near
#the actual saturation temperature, so we use the direct algorithm. otherwise, we use a safe approach, starting from the critical
#coordinate and descending.

#by default, starts right before the critical point, and descends via Clapeyron equation: (∂p/∂T)sat = ΔS/ΔV ≈ Δp/ΔT

function saturation_temperature_impl(model::EoSModel,p,method::ClapeyronSaturation)
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

export ClapeyronSaturation