"""
    x0_volume_liquid(model,T,z)

Returns an initial guess to the liquid volume, dependent on temperature and composition. by default is 1.25 times [`lb_volume`](@ref).
"""
function x0_volume_liquid(model,T,z)
    v_lb = lb_volume(model,z)
    return v_lb*1.25
end

"""
    x0_volume_gas(model,p,T,z)

Returns an initial guess to the gas volume, depending of pressure, temperature and composition. by default uses [`volume_virial`](@ref)
"""
function x0_volume_gas(model,p,T,z)
    return volume_virial(model,p,T,z)
end

"""
    x0_volume_solid(model,T,z)

Returns an initial guess to the solid volume, dependent on temperature and composition. needs to be defined for EoS that support solid phase. by default returns NaN
"""
function x0_volume_solid(model,T,z)
    _0 = zero(T+first(z))
    return _0/_0
end

"""
    x0_volume(model,p,T,z; phase = :unknown)

Returns an initial guess of the volume at a pressure, temperature, composition and suggested phase.

If the suggested phase is `:unkwown` or `:liquid`, calls [`x0_volume_liquid`](@ref).

If the suggested phase is `:gas`, calls [`x0_volume_gas`](@ref).

If the suggested phase is `solid`, calls [`x0_volume_solid`](@ref).

Returns `NaN` otherwise

"""
function x0_volume(model, p, T, z = SA[1.0]; phase = :unknown)
    phase = Symbol(phase)
    if phase === :unknown || is_liquid(phase)
        return x0_volume_liquid(model,T,z)
    elseif is_vapour(phase)
        return x0_volume_gas(model,p,T,z)
    elseif is_supercritical(phase)
        return x0_volume_gas(model,p,T,z)
    elseif is_solid(phase)
        x0_volume_solid(model,T,z)
    else
        _0 = zero(p+T+first(z))
        return _0/_0
    end
end

"""
    lb_volume(model::EoSModel,z=SA[1.0])

Returns the lower bound volume.

It has different meanings depending on the Equation of State, but symbolizes the minimum allowable volume at a certain composition:

- SAFT EoS: the packing volume
- Cubic EoS, covolume (b) parameter

On empiric equations of state, the value is chosen to match the volume of the conditions at maximum pressure and minimum temperature
, but the equation itself normally can be evaluated at lower volumes.

On SAFT and Cubic EoS, volumes lower than `lb_volume` will likely error.

The lower bound volume is used for guesses of liquid volumes at a certain pressure, saturated liquid volumes and critical volumes.
"""
function lb_volume end

"""
    T_scale(model::EoS,z=SA[1.0])

Represents a temperature scaling factor.

On any EoS based on Critical parameters (Cubic or Empiric EoS), the temperature scaling factor is chosen to be the critical temperature.

On SAFT or other molecular EoS, the temperature scaling factor is chosen to be a function of the potential depth ϵ.

Used as scaling factors in [`saturation_pressure`](@ref) and as input for solving [`crit_pure`](@ref)
"""
function T_scale end

"""
    p_scale(model::SAFTModel,z=SA[1.0])

Represents a pressure scaling factor

On any EoS based on Critical parameters (Cubic or
Empiric EoS), the pressure scaling factor is
chosen to be a function of the critical pressure.

On SAFT or other molecular EoS, the temperature
scaling factor is chosen to a function of ∑(zᵢ*ϵᵢ*(σᵢᵢ)³)

Used as scaling factors in [`saturation_pressure`](@ref) and as input for solving [`crit_pure`](@ref)

"""
function p_scale end

"""
    antoine_coef(model)

should return a 3-Tuple containing reduced Antoine Coefficients. The Coefficients follow the correlation:
```
lnp̄ = log(p / p_scale(model))
T̃ = T/T_scale(model)
lnp̄ = A - B/(T̄ + C))
```

By default returns `nothing`. This is to use alternative methods in case Antoine coefficients aren't available. Used mainly in single and multicomponent temperature calculations.
"""
function antoine_coef end

antoine_coef(model) = nothing


"""
    x0_sat_pure(model::EoSModel,T,z=SA[1.0])

Returns a 2-tuple corresponding to `(Vₗ,Vᵥ)`, where `Vₗ` and `Vᵥ` are the liquid and vapor initial guesses.

Used in [`saturation_pressure`](@ref) methods that require initial volume guesses.

It can be overloaded to provide more accurate estimates if necessary.
"""
function x0_sat_pure(model,T,z=SA[1.0])
    #=theory as follows
    #given T = Teos:
    #calculate B(Teos)
    PB = Peos = -2B

    veos = volume_compress(model,PB,T)
    with a (P,V,T,B) pair, change to
    (P,V,a,b)
    =#


    B = second_virial_coefficient(model,T,SA[1.0])
    lb_v = lb_volume(model,SA[1.0])*one(T)
    _0 = zero(B)
    #virial volume below lower bound volume.
    #that means that we are way over the critical point
    if -2B < lb_v
        _nan = _0/_0
        return (_nan,_nan)
    end

    p = -0.25*R̄*T/B
    vl = x0_volume(model,p,T,z,phase=:l)
    vl = _volume_compress(model,p,T,SA[1.0],vl)
    #=the basis is that p = RT/v-b - a/v2
    we have a (p,v,T) pair
    and B = 2nd virial coefficient = b-a/RT
    with that, we solve for a and b
    as a and b are vdW, Pc and Tc can be calculated
    with Tc and Pc, we use [1] to calculate vl0 and vv0
    with Tc, we can also know in what regime we are.
    in near critical pressures, we use directly vv0 = -2B
    and vl0 = 4*lb_v
    [1]
    DOI: 10.1007/s10910-007-9272-4
    Journal of Mathematical Chemistry, Vol. 43, No. 4, May 2008 (© 2007)
    The van der Waals equation: analytical and approximate solutions
    =#
    γ = p*vl*vl/(R̄*T)
    _c = vl*(vl + B - γ)
    _b = γ - B - vl
    Δ = Solvers.det_22(_b,_b,4,_c)
    if isnan(vl) | (Δ < 0)
        #fails on two ocassions:
        #near critical point, or too low.
        #old strategy
        x0l = 4*lb_v
        x0v = -2*B + 2*lb_v
        return (x0l,x0v)
    end
    Δsqrt = sqrt(Δ)
    b1 = 0.5*(-_b + Δsqrt)
    b2 = 0.5*(-_b - Δsqrt)
    if b1 < 0
        b = b2
    elseif b1 > vl
        b = b2
    else
        b = b1
    end
    a = -R̄*T*(B-b)
    Vc = 3*b
    Ωa =  27/64
    Ωb =  1/8
    ar = a/Ωa
    br = b/Ωb
    Tc = ar/br/R̄
    Pc = ar/(br*br)
    Tr = T/Tc
    #Tr(vdW approx) > Tr(model), try default.
    if Tr >= 1
        x0l = 4*lb_v
        x0v = -2*B + 2*lb_v
        return (x0l,x0v)
    end
    #@show (Tc,Pc)
    # if b1/b2 < -0.95, then T is near Tc.
    #if b<lb_v then we are in trouble
    #critical regime or something horribly wrong happened
    if (b1/b2 < -0.95) | (b<lb_v) | (Tr>0.99)
        println("warning: near criticality")
        x0l = 4*lb_v
        x0v = -2*B #gas volume as high as possible
        return (x0l,x0v)
    end
    Vl0,Vv0 = vdw_x0_xat_pure(T,Tc,Pc,Vc)
    x0l = min(Vl0,vl)
    x0v = min(1e4*one(Vv0),Vv0) #cutoff volume
    return (x0l,x0v)
end

function vdw_x0_xat_pure(T,T_c,P_c,V_c)
    Tr = T/T_c
    Trm1 = 1.0-Tr
    Trmid = sqrt(Trm1)
    if Tr >= 0.7
        c_l = 1.0+2.0*Trmid + 0.4*Trm1 - 0.52*Trmid*Trm1 +0.115*Trm1*Trm1 #Eq. 29
    else
        c_l = 1.5*(1+sqrt(1-(32/27)*Tr)) #Eq. 32
    end

    if Tr >= 0.46
        #Eq. 30, valid in 0.46 < Tr < 1
        c_v = 1.0-2.0*Trmid + 0.4*Trm1 + 0.52*Trmid*Trm1 +0.207*Trm1*Trm1
    elseif Tr <= 0.33
        #Eq. 33, valid in 0 < Tr < 0.33
        c_v = (3*c_l/(ℯ*(3-c_l)))*exp(-(1/(1-c_l/3)))

    else
        #Eq. 31 valid in 0.25 < Tr < 1
        mean_c = 1.0 + 0.4*Trm1 + 0.161*Trm1*Trm1
        c_v = 2*mean_c - c_l
    end
    #volumes predicted by vdW
    Vl0 = V_c/c_l
    Vv0 = V_c/c_v
    return (Vl0,Vv0)
end

function scale_sat_pure(model,z=SA[1.0])
    p    = 1/p_scale(model,z)
    μ    = 1/R̄/T_scale(model,z)
    return p,μ
end

"""
    x0_psat(model::EoSModel, T,crit = nothing)

Initial point for saturation pressure, given the temperature and V,T critical coordinates.
On moderate pressures it will use a Zero Pressure initialization. On pressures near the critical point it will switch to spinodal finding.

Used in [`saturation_pressure`](@ref) methods that require initial pressure guesses.

if the initial temperature is over the critical point, it returns `NaN`.

It can be overloaded to provide more accurate estimates if necessary.
"""
function x0_psat(model,T,crit = nothing)
    coeffs = antoine_coef(model)
    if coeffs !== nothing
        A,B,C = coeffs
        T̃ = T/T_scale(model)
        lnp̃ = A - B/(T̃ + C)
        ps = p_scale(model)
        px =  exp(lnp̃)*ps
        return px
    end
    if isnothing(crit)
        crit = crit_pure(model)
    end
    Tc, Pc, Vc = crit
    if T > Tc
        return zero(T)/zero(T)
    end
    return x0_psat(model, T, Tc, Vc)
end

function x0_psat(model::EoSModel, T, Tc, Vc)
    # Function to get an initial guess for the saturation pressure at a given temperature
    z = SA[1.] #static vector
    _0 = zero(T+Tc+Vc)
    RT = R̄*T
    Tr = T/Tc
    # Zero pressure initiation
    if Tr < 0.8
        P0 = _0
        vol_liq0 = volume(model, P0, T, phase=:liquid)
        ares = a_res(model, vol_liq0, T, z)
        lnϕ_liq0 = ares - 1. + log(RT/vol_liq0)
        P0 = exp(lnϕ_liq0)
    # Pmin, Pmax initiation
    elseif Tr <= 1.0
        low_v = Vc
        up_v = 5 * Vc
        #note: P_max is the pressure at the maximum volume, not the maximum pressure
        fmax(V) = -pressure(model, V, T)
        sol_max = Solvers.optimize(fmax, (low_v, up_v))
        P_max = -Solvers.x_minimum(sol_max)
        low_v = lb_volume(model)
        up_v = Vc
        #note: P_min is the pressure at the minimum volume, not the minimum pressure
        fmin(V) = pressure(model, V, T)
        sol_min = Solvers.optimize(fmin, (low_v,up_v))
        P_min = Solvers.x_minimum(sol_min)
        P0 = (max(zero(P_min), P_min) + P_max) / 2
    else
        P0 = _0/_0 #NaN, but propagates the type
    end
    return P0
end

"""
    x0_saturation_temperature(model::EoSModel,p)

Returns a 3-tuple corresponding to `(T,Vₗ,Vᵥ)`, `T` is the initial guess for temperature and `Vₗ` and `Vᵥ` are the liquid and vapor initial guesses.
Used in [`saturation_temperature`](@ref) with [`AntoineSaturation`](@ref).
"""
function x0_saturation_temperature end

function x0_saturation_temperature(model::EoSModel,p,::Nothing)
    crit = crit_pure(model)
    return x0_saturation_temperature(model,p,crit)
end

function x0_saturation_temperature(model::EoSModel,p,crit::Tuple)
    Tc,Pc,Vc = crit
    A,B,C = (6.668322465137264,6.098791871032391,-0.08318016317721941)
    if Pc < p
        nan = zero(p)/zero(p)
        return (nan,nan,nan)
    end
    lnp̄ = log(p / Pc)
    T0 = Tc*(B/(A-lnp̄)-C)
    pii,vli,vvi = saturation_pressure(model,T0,ChemPotVSaturation(;crit))

    if isnan(pii)
        nan = zero(p)/zero(p)
        return (nan,nan,nan)
    end

    Δp = (p-pii)
    S_v = VT_entropy(model,vvi,T0)
    S_l = VT_entropy(model,vli,T0)
    ΔS = S_v - S_l
    ΔV = vvi - vli
    dpdt = ΔS/ΔV #≈ (p - pii)/(T-Tnew)
    T = T0 + Δp/dpdt
    vv = volume_virial(model,p,T)
    vl = 0.3*lb_volume(model) + 0.7*vli
    return (T,vl,vv)
end

"""
    x0_crit_pure(model::EoSModel)

Returns a 2-tuple corresponding to
    `(k,log10(Vc0))`, where `k` is `Tc0/T_scale(model,z)`
"""
function x0_crit_pure end

function x0_crit_pure(model::EoSModel)
    lb_v = lb_volume(model)
    (1.5, log10(lb_v/0.3))
end

#=
the following methods are fallbacks,
that require just the definition of T_scale,p_scale and lb_volume
respectively. if possible, each eos should define those
=#
function T_scales(model)
    n = length(model)
    x = zeros(n)
    res = zeros(n)
    for i = 1:n
        x[i] = 1.0
        res[i] = T_scale(model,x)
        x[i] = 0.0
    end
    return res
end

T_scales(model,z) = T_scales(model)

function x0_sat_pure2(model,T,z=SA[1.0])
    #=
    x0_sat_pure v3: full local van der waals.

    first version: 
        - constant values
    second version: 
        - calculate B(Teos) (B = b - a/RT)
        - PB = Peos = P(v = -2B,T) # gas spinodal aproximation
        - vl obtained from volume(model,PB,T)
        - solve a,b from (P,V,T,B) coordinates 
        - in case that volume solver fails, regress to constant near critical.
    third version:
        - calculate B(Teos) (B = b - a/RT)
        - solve (a,b) such as:
            - p(eos) = max(p(vdw),0)
            - v ≈ b(1 + sqrt(0.5RTb/a)) #liquid spinodal aproximation
            - (modified volume_compress, does not fail even near criticality.)
    =#
    B = second_virial_coefficient(model,T,SA[1.0])
    lb_v = lb_volume(model,SA[1.0])*one(T)
    _0 = zero(B)
    TT = typeof(_0)
    if -2B < lb_v 
        _nan = _0/_0
        return (_nan,_nan)
    end
    
    a,b = vdw_ab(model,T,false,B)
    #@show a,b
    Tc,Pc,Vc = crit_pure_vdw_ab(a,b)

    #=the basis is that p = RT/v-b - a/v2
    we have a (p,v,T) pair
    and B = 2nd virial coefficient = b-a/RT
    with that, we solve for a and b 
    as a and b are vdW, Pc and Tc can be calculated
    with Tc and Pc, we use [1] to calculate vl0 and vv0
    
    with Tc, we can also know in what regime we are.
    in near critical pressures, we use directly vv0 = -2B
    and vl0 = 4*lb_v

    [1]
    DOI: 10.1007/s10910-007-9272-4
    Journal of Mathematical Chemistry, Vol. 43, No. 4, May 2008 (© 2007)
    The van der Waals equation: analytical and approximate solutions
    =#

    #pseudo spinodal values
   
    v_sv = -2*B
    p_sv = -0.25*R̄*T/B
    v_sl = b*(1+sqrt(0.5*R̄*T*b/a))
    p_sl =  max(R̄*T/(v_sl- b) - a/(v_sl^2),zero(v_sl))

    Tr = T/Tc
    #Tr(vdW approx) > Tr(model)
    if Tr >= 1
        _nan = _0/_0
        return (_nan,_nan)
    end

    if p_sv < p_sl ||  p_sl > Pc 
        #if p_sl is over the critical point or over the vapour spinodal
        #then we are in a critical regime. we then recalculate a,b 
        #with a critical correlation. 
        a,b = vdw_ab(model,T,true,B,b)
        Tc,Pc,Vc = crit_pure_vdw_ab(a,b)
        v_sv = -2*B
        p_sv = -0.25*R̄*T/B
        v_sl = b*(1+sqrt(0.5*R̄*T*b/a))
        p_sl =  max(R̄*T/(v_sl- b) - a/(v_sl^2),zero(v_sl))
        #use rackett with the critical volume
        Zc = Pc*Vc/(R̄*Tc)
        vl = (R̄*Tc/Pc)*Zc^(1 + (1-Tr)^(2/7))
        vv = v_sv
        vlc,vvc =  vdw_x0_xat_pure(T,Tc,Pc,Vc)
        vl = 0.5*(vl+vlc)
        vv = 0.5*(vv+vvc)
    else
        p_l = p_sl
        pm = 0.5*(p_sv+p_l)
        #one volume iteration
        dpdv_sl = -R̄*T/(v_sl - b)^2 + 2*a/v_sl^3
        _Δ = (pm-p_l)/(v_sl*dpdv_sl)
        vl = v_sl*exp(_Δ)
        vv = volume_virial(B,pm,T)
    end
    return (log10(vl),log10(vv)) 
end

function vdw_ab(model,T,crit = false,B = second_virial_coefficient(model,T),b0 = x0_volume_liquid(model,T,SA[1.0]))
    z = SA[1.0]
    lb_v = lb_volume(model,z)
    x0_l = x0_volume_liquid(model,T,z)
    logb0 = log(x0_l)
    log_lb_v = log(lb_v)
    _nan = zero(T)/zero(T)
    _0 = zero(_nan)
    Ωa =  27/64
    Ωb =  1/8
    function logstep(logb)       
        logb < log_lb_v && return _nan
        b = exp(logb)
        a = -R̄*T*(B - b)
        if crit
            Tc,Pc,Vc = crit_pure_vdw_ab(a,b)
            Tr = T/Tc
            Trm1 = 1.0-Tr
            Trmid = sqrt(Trm1)
            c_l = 1.0 + 2.0*Trmid + 0.4*Trm1 - 0.52*Trmid*Trm1 + 0.115*Trm1*Trm1
            v = 0.5*(b*(1+sqrt(0.5*R̄*T*b/a)) + (Vc/c_l))
        else
            v = b*(1+sqrt(0.5*R̄*T*b/a))
        end
        p_vdw = max(R̄*T/(v- b) - a/(v^2),_0)
        p_eos,dpdV = p∂p∂V(model,v,T,z)
        dpdV > 0 && return _nan #inline mechanical stability. 
        abs(p_eos-p_vdw) < 3eps(p_eos) && return _0
        Δ = (p_vdw-p_eos)/(v*dpdV) 
        return Δ
    end

    function f_fixpoint(_V)
        Δ = logstep(_V)
        vv = _V + Δ
        return vv
    end
    log_bb = @nan(Solvers.fixpoint(f_fixpoint,logb0,Solvers.SSFixPoint(),rtol = 1e-8,max_iters=100),_nan)
    bb =  exp(log_bb)
    aa = -R̄*T*(B - bb)
    return (aa,bb)
end

function crit_pure_vdw_ab(a,b)
    Ωa =  27/64
    Ωb =  1/8
    ar = a/Ωa
    br = b/Ωb
    Tc = ar/br/R̄
    Pc = ar/(br*br)
    Vc = 3*b
    return Tc,Pc,Vc
end
