"""
    x0_volume_liquid(model,T,z)

Returns an initial guess to the liquid volume, dependent on temperature and composition. by default is 1.25 times the [lower bound volume](@ref lb_volume).
"""
function x0_volume_liquid(model,T,z)
    v_lb = lb_volume(model,z)
    return v_lb*1.25
end

"""
    x0_volume_gas(model,p,T,z)

Returns an initial guess to the gas volume, depending of pressure, temperature and composition. by default uses a [virial aproximation](@ref volume_virial)
"""
function x0_volume_gas(model,p,T,z)
    return volume_virial(model,p,T,z)
end

"""
    x0_volume(model,p,T,z; phase = :unknown)

Returns an initial guess of the volume at a pressure, temperature, composition and suggested phase.

If the suggested phase is `:unkwown` or `:liquid`, calls [`x0_volume_liquid`](@ref).

If the suggested phase is `:gas`, calls [`x0_volume_gas`](@ref).

"""
function x0_volume(model,p,T,z; phase = :unknown)
    phase = Symbol(phase)
    if phase === :unknown || is_liquid(phase)
        return x0_volume_liquid(model,T,z)
    elseif is_vapour(phase)
        return x0_volume_gas(model,p,T,z)
    elseif is_supercritical(phase)
     else
        error("unreachable state on x0_volume")
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
    x0_sat_pure(model::EoSModel,T,z=SA[1.0])

Returns a 2-tuple corresponding to `(log10(Vₗ),log10(Vᵥ))`, where Vₗ and Vᵥ are the liquid and vapor initial guesses. 
Used in [`sat_pure`](@ref).
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
    _c = vl*vl + B*vl - γ*vl 
    _b = γ - B - vl
    Δ = _b*_b - 4*_c
    if isnan(vl) | (Δ < 0)
        #fails on two ocassions:
        #near critical point, or too low.
        #old strategy
        x0l = 4*lb_v
        x0v = -2*B + 2*lb_v
        return (log10(x0l),log10(x0v))
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
    #Tr(vdW approx) > Tr(model)
    if Tr >= 1
        _nan = _0/_0
        return (_nan,_nan)
    end
    # if b1/b2 < -0.95, then T is near Tc.
    #if b<lb_v then we are in trouble 
    #critical regime or something horribly wrong happened
    if (b1/b2 < -0.95) | (b<lb_v) | (Tr>0.99)
        x0l = 4*lb_v
        x0v = -2*B #gas volume as high as possible
        return (log10(x0l),log10(x0v))   
    end
    Vl0,Vv0 = vdw_x0_xat_pure(T,Tc,Pc,Vc)
    x0l = min(Vl0,vl)
    x0v = min(1e4*one(Vv0),Vv0) #cutoff volume
    return (log10(x0l),log10(x0v)) 
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
    Vl0 = (1/c_l)*V_c
    Vv0 = (1/c_v)*V_c
    return (Vl0,Vv0)
end

function scale_sat_pure(model,z=SA[1.0])
    p    = 1/p_scale(model,z)
    μ    = 1/R̄/T_scale(model,z)
    return p,μ
end

"""
    x0_crit_pure(model::SAFTModel)

Returns a 2-tuple corresponding to
    `(k,log10(Vc0))`, where `k` is `Tc0/T_scale(model,z)`
"""
function x0_crit_pure end

function x0_crit_pure(model::EoSModel)
    lb_v = lb_volume(model)
    (1.5, log10(lb_v/0.3))
end

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
