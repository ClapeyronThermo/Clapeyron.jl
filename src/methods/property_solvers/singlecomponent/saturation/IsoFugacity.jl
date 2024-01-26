"""
    IsoFugacitySaturation <: SaturationMethod
    IsoFugacitySaturation(;p0 = nothing,
        vl = nothing,
        vv = nothing,
        crit = nothing,
        max_iters = 20,
        p_tol = sqrt(eps(Float64)))

Saturation method for `saturation_pressure`. Uses the isofugacity criteria. Ideal for Cubics or other EoS where the volume calculations are cheap. 
If `p0` is not provided, it will be calculated via [`x0_psat`](@ref).
"""
struct IsoFugacitySaturation{T,V,C} <: SaturationMethod
    p0::T
    vl::V
    vv::V
    crit::C
    max_iters::Int
    p_tol::Float64
end

function IsoFugacitySaturation(;p0 = nothing,
                                vl = nothing,
                                vv = nothing,
                                crit = nothing,
                                max_iters = 20,
                                p_tol = sqrt(eps(Float64)))
    p0 === nothing && (p0 = NaN)
    
    if (vv !== nothing) & (vl !== nothing)
        p0,vl,vv = promote(p0,vl,vv)
        return IsoFugacitySaturation(p0,vl,vv,crit,max_iters,p_tol)
    elseif vv === nothing && vl === nothing
        return IsoFugacitySaturation(p0,vl,vv,crit,max_iters,p_tol)
    else
        throw(error("you need to specify both vl and vv"))
    end
end

function saturation_pressure_impl(model::EoSModel,T,method::IsoFugacitySaturation)
    vol0 = (method.vl,method.vv,T)
    p0 = method.p0
    if isnan(p0)
        p0 = x0_psat(model, T, method.crit)
    end

    if isnan(p0) #over critical point, or something else.
        nan = p0/p0
        return (nan,nan,nan)
    end

    return psat_fugacity(model,T,p0,vol0,method.max_iters,method.p_tol)
end

function psat_fugacity(model::EoSModel, T, p0, vol0=(nothing, nothing),max_iters = 20,p_tol = sqrt(eps(Float64)))
    # Objetive function to solve saturation pressure using the pressure as iterable variable
    # T = Saturation Temperature
    # p0 = initial guess for the saturation pressure
    # vol0 = initial guesses for the phase volumes = [vol liquid, vol vapor]
    # out = Saturation Pressure, vol liquid, vol vapor
    z = SA[1.]
    RT = Rgas(model)*T
    P = 1. * p0
    vol_liq0, vol_vap0 = vol0
    #we use volume here, because cubics can opt in to their root solver.
    vol_liq0 === nothing && (vol_liq0 = volume(model,P,T,z,phase =:liquid))
    vol_vap0 === nothing && (vol_vap0 = volume(model,P,T,z,phase =:gas))
    
    vol_liq = vol_liq0 
    vol_vap = vol_vap0
    #@show vol_liq, vol_vap
    itmax = max_iters
    fun(_V) = eos_res(model, _V, T,SA[1.])

    for i in 1:itmax
        # Computing chemical potential
        A_l,Av_l = Solvers.f∂f(fun,vol_liq)
        A_v,Av_v = Solvers.f∂f(fun,vol_vap)
        μ_liq = muladd(-vol_liq,Av_l,A_l)
        μ_vap = muladd(-vol_vap,Av_v,A_v)
        #μ_liq = VT_chemical_potential_res(model, vol_liq, T)[1]
        #μ_vap = VT_chemical_potential_res(model, vol_vap, T)[1]

        Z_liq = P*vol_liq/RT
        Z_vap = P*vol_vap/RT
        
        lnϕ_liq = μ_liq/RT - log(Z_liq)
        lnϕ_vap = μ_vap/RT - log(Z_vap)
        # Updating the saturation pressure
        FO = lnϕ_vap - lnϕ_liq
        dFO = (Z_vap - Z_liq) / P
        #dP = -expm1(lnϕ_liq-lnϕ_vap)
        dP = FO / dFO
        #at very low pressures and some EoS (GERG2008(["water"]) at 298K) the first step could be too big.
        #this will limit steps diving on negative pressures
        P = max(P - dP,P/2)
        
        if abs(dP) < p_tol; break; end
        # Updating the phase volumes
        vol_liq = _volume_compress(model, P, T, z,vol_liq)
        vol_vap = _volume_compress(model, P, T, z,vol_vap)
    end
    return P, vol_liq, vol_vap
end

export IsoFugacitySaturation
