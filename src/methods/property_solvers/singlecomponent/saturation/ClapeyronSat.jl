struct ClapeyronSaturation{T} <: SaturationMethod
    Temp::Union{Nothing,T}
end

ClapeyronSaturation() = ClapeyronSaturation{Nothing}(nothing)
#if a model overloads x0_saturation_temperature to return a T0::Number, we can assume this number is near
#the actual saturation temperature, so we use the direct algorithm. otherwise, we use a safe approach, starting from the critical
#coordinate and descending.

#by default, starts right before the critical point, and descends via Clapeyron equation: (∂p/∂T)sat = ΔS/ΔV ≈ Δp/ΔT

function saturation_temperature_impl(model::EoSModel,p,method::ClapeyronSaturation{Nothing})
    Tc,pc,vc = crit_pure(model)
    TT = typeof(vc*pc/Tc)
    nan = zero(TT)/zero(TT)
    p > 0.99999pc && (return (nan,nan,nan)) 
    T0 = 0.99*Tc
    isnan(T0) && (return (nan,nan,nan))
    method_init = ClapeyronSaturation(T0)
    return saturation_temperature_impl(model,p,method_init)
end

function saturation_temperature_impl(model::EoSModel,p,method::ClapeyronSaturation)
    T0 = method.Temp/one(method.Temp)
    TT = typeof(T0)
    nan = zero(T0)/zero(T0)
    cache = Ref{Tuple{TT,TT,TT,TT,Bool}}((nan,nan,nan,nan,false))
    f(T) = Obj_sat_pure_T(model,T,p,cache)
    T = Solvers.fixpoint(f,T0)
    _,_,v_l,v_v,_ = cache[]
    return T,v_l,v_v
end

function Obj_sat_pure_T(model,T,p,cache)
    Told,pold,vlold,vvold,use_v = cache[]
    if use_v  
        sat_method  = ChemPotVSaturation(log10(vlold*0.99),log10(1.01*vvold))
    else
        sat_method = ChemPotVSaturation(x0_sat_pure(model,T))
    end
    pii,vli,vvi = saturation_pressure(model,T,sat_method)
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
    cache[] = (T,pii,vli,vvi,use_v) 
    S_v = VT_entropy(model,vvi,T)
    S_l = VT_entropy(model,vli,T)
    ΔS = S_v - S_l
    ΔV = vvi - vli
    dpdt = ΔS/ΔV #≈ (p - pii)/(T-Tnew)
    Ti = T + Δp/dpdt
    return Ti
end
#=
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
=#
export ClapeyronSaturation