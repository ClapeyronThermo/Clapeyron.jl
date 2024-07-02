"""
    x0_volume_liquid(model,T,z)
    x0_volume_liquid(model,p,T,z)

Returns an initial guess to the liquid volume, dependent on temperature and composition. by default is 1.25 times [`lb_volume`](@ref).
"""
function x0_volume_liquid(model,T,z)
    v_lb = lb_volume(model,z)
    return v_lb*1.25
end

x0_volume_liquid(model,p,T,z) = x0_volume_liquid(model,T,z)

"""
    x0_volume_gas(model,p,T,z)

Returns an initial guess to the gas volume, depending of pressure, temperature and composition. by default uses [`volume_virial`](@ref)
"""
function x0_volume_gas(model,p,T,z)
    return volume_virial(model,p,T,z)
end

"""
    x0_volume_solid(model,T,z)
    x0_volume_solid(model,p,T,z)

Returns an initial guess to the solid volume, dependent on temperature and composition. needs to be defined for EoS that support solid phase. by default returns NaN. can be overrided if the EoS defines `is_solid(::EoSModel) = true`
"""
function x0_volume_solid(model,T,z)
    if is_solid(model)
        v_lb = lb_volume(model,z)
        return v_lb*1.05
    else
        _0 = zero(T+first(z))
        return _0/_0
    end
end

x0_volume_solid(model,p,T,z) = x0_volume_solid(model,T,z)

"""
    x0_volume(model,p,T,z; phase = :unknown)
Returns an initial guess of the volume at a pressure, temperature, composition and suggested phase.
If the suggested phase is `:unknown` or `:liquid`, calls [`x0_volume_liquid`](@ref).
If the suggested phase is `:gas`, calls [`x0_volume_gas`](@ref).
If the suggested phase is `solid`, calls [`x0_volume_solid`](@ref).
Returns `NaN` otherwise
"""
function x0_volume(model, p, T, z = SA[1.0]; phase = :unknown)
    return x0_volume_impl(model,p,T,z,phase)
end

function x0_volume_impl(model, p, T, z = SA[1.0], phase = :unknown)
    phase = Symbol(phase)
    if is_unknown(equilibrium) || is_liquid(phase)
        return x0_volume_liquid(model,p,T,z)
    elseif is_vapour(phase)
        return x0_volume_gas(model,p,T,z)
    elseif is_supercritical(phase)
        return x0_volume_gas(model,p,T,z)
    elseif is_solid(phase)
        x0_volume_solid(model,p,T,z)
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
    x0_sat_pure(model::EoSModel,T)
Returns a 2-tuple corresponding to `(Vₗ,Vᵥ)`, where `Vₗ` and `Vᵥ` are the liquid and vapor initial guesses.
Used in [`saturation_pressure`](@ref) methods that require initial volume guesses.
It can be overloaded to provide more accurate estimates if necessary.
"""
x0_sat_pure(model,T) = x0_sat_pure_virial(model,T)

function x0_sat_pure_virial(model,T)
    z = SA[1.0]
    single_component_check(x0_sat_pure,model)

    #=theory as follows
    #given T = Teos:
    #calculate B(Teos)
    PB = Peos = -2B
    veos = volume_compress(model,PB,T)
    with a (P,V,T,B) pair, change to
    (P,V,a,b)
    =#

    R̄ = Rgas(model)
    RT= R̄*T
    B = second_virial_coefficient(model,T,z)
    _0,_1 = zero(B),oneunit(B)
    lb_v = lb_volume(model,z)*_1
    #=
    some very complicated models, like DAPT, fail on the calculation of the second virial coefficient.
    while this is a numerical problem in the model itself, it is better to catch this early.

    while a volume calculation is harder
    =#
    if isnan(B)
        x0l = 3*lb_v
        px = pressure(model,x0l,T)
        if px < 0 #low pressure
            return x0_sat_pure_near0(model,T)
        else #high pressure?
            return 4*lb_v,20*lb_v
        end
    end


    _0 = zero(B)
    #virial volume below lower bound volume.
    #that means that we are way over the critical point
    if -2B < lb_v
        _nan = _0/_0
        return (_nan,_nan)
    end

    #=
    we define the following functions
    γc = Pc(eos)/P(virial, at v = -2B, B = B(eos,Tc))
    γT(T) = P(v = -2B,T)/P(virial, at v = -2B, B = B(eos,T))

    for cubics, γc is constant
    we can relate γc = γc(γT(Tr = 1)), to obtain a pressure that has a liquid root.
    we can further extend this to γc = γc(γT(T))
    because at near critical pressures, the virial predicted pressure is below the liquid spinodal pressure
    in one sense, γc is a correction factor.
    =#
    pv_virial = -0.25*RT/B #maximum virial predicted pressure
    vv_virial = -2*B #volume at maximum virial pressure
    pv_eos = pressure(model,vv_virial,T) #eos predicted pressure, gas phase

    γT = pv_eos/pv_virial

    #fitted function, using all coolprop fluids, at Tr = 1
    aγ,bγ,cγ = 1.2442071971165476e-5, -8.695786307570637, 1.0505452946870144
    γc = aγ*exp(-γT*bγ) + cγ
    pl0 = γc*pv_virial #pressure of which we are (almost) sure there exists a liquid root
    vl_x0 = x0_volume(model,pl0,T,z,phase=:l)
    vl = _volume_compress(model,pl0,T,z,vl_x0)
    #=the basis is that p = RT/v-b - a/v2
    we can interpolate a vdW EoS between a liquid and a gas (p,v) point.
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
    a,b = vdw_coefficients(vl,pl0,vv_virial,pv_eos,T)
    Tc,Pc,Vc = vdw_crit_pure(a,b)
    if isnan(a) | isnan(b)
        #fails on two ocassions:
        #near critical point, or too low.
        #return high pressure estimate
        return 4*lb_v,vv_virial + 2*lb_v
    end
    Tr = T/Tc
    #Tr(vdW approx) > Tr(model), try default.
    if Tr >= 1
        x0l = 4*lb_v
        x0v = vv_virial + 2*lb_v
        return (x0l,x0v)
    end
    #saturation volumes calculated from a vdW EoS.
    Vl0,Vv0 = vdw_x0_sat_pure(T,Tc,Pc,Vc)

    if Tr > 0.97
        #in this region, finding both spinodals is cheaper than incurring in a crit_pure calculation.
        return x0_sat_pure_hermite_spinodal(model,Vl0,Vv0,T)
    elseif Vv0 > 1e4*one(Vv0)
        #gas volume over threshold. but not diverged.
        #normally this happens at low temperatures. we could suppose that Vl0 is a
        #"zero-pressure" volume, apply corresponding strategy
        x0l = min(Vl0,vl)
        return x0_sat_pure_near0(model,T,x0l)
    else
        #we correct the gas saturation pressure. normally, B_vdw is lower than B_eos.
        #this causes underprediction on the vapour initial point.
        B_vdw = b - a/RT
        p_vdw_sat = RT/(Vv0 - b) - a/Vv0/Vv0
        p_vdw_gas_corrected = p_vdw_sat*B/B_vdw
        pv_vdw_virial = -0.25*RT/B_vdw
        Vv02 = volume_virial(B_vdw,p_vdw_gas_corrected,T) # virial correction.
        if pv_virial < p_vdw_gas_corrected
            return x0_sat_pure_hermite_spinodal(model,Vl0,Vv02,T)
        end
        return Vl0,Vv02
    end
end

function vdw_coefficients(vl,pl,vv,pv,T)
    #px = RT/(vx-b) - a/vx/vx
    RT = R̄*T
    βl = pl*vl*vl/RT # vl*vl/(vl - b) - a/RT
    βv = pv*vv*vv/RT # vv*vv/(vv - b) - a/RT
    Δβ = βv - βl
    vv2,vl2 = vv*vv,vl*vl
    vvΔβ,vlΔβ = vv2/Δβ,vl2/Δβ
    #=
    #we solve a in terms of the liquid values
    a/RT = vl*vl/(vl - b) - βl

    then:

    βv - βl = vv2/(vv - b) - vl2/(vl - b)
    Δβ*(vl - b)*(vv - b) = vv2*(vl - b) - vl2*(vv - b)
    (vl - b)*(vv - b) = vv2/Δβ*(vl - b) - vl2/Δβ*(vv - b)
    (vl - b)*(vv - b) = vvΔβ*(vl - b) - vlΔβ*(vv - b)
    b^2 -(vl + vv)*b + vl*vv = vvΔβ*vl - vvΔβ*b - vlΔβ*vv + vlΔβ*b
    b^2 -(vl + vv)*b + vl*vv - vvΔβ*vl + vvΔβ*b + vlΔβ*vv - vlΔβ*b = 0
    b^2 + (-vl + -vv + vvΔβ - vlΔβ)*b + (vl*vv - vvβ*vl + vlΔβ*vv) = 0
    =#
    _c = vl*vv - vvΔβ*vl + vlΔβ*vv
    _b = -vl + -vv + vvΔβ - vlΔβ

    #quadratic solver
    Δ = Solvers.det_22(_b,_b,4,_c)
    if isnan(vl) | (Δ < 0)
        nan = zero(Δ)/zero(Δ)
        return nan,nan
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
    a = RT*(vl2/(vl - b) - βl)
    return a,b
end

#vdw crit pure values,given a and b
function vdw_crit_pure(a,b)
    Ωa,Ωb = 27/64, 1/8
    Vc = 3*b
    ar = a/Ωa
    br = b/Ωb
    Tc = ar/br/R̄
    Pc = ar/(br*br)
    return Tc,Pc,Vc
end

function x0_sat_pure_near0(model, T,vl0 = volume(model,zero(T),T,phase =:liquid))
    R̄ = Rgas(model)
    z = SA[1.0]
    RT = R̄*T
    ares = a_res(model, vl0, T, z)
    lnϕ_liq0 = ares - 1 + log(RT/vl0)
    P = exp(lnϕ_liq0)
    vl = volume(model,P,T,z,vol0 = vl0)
    vv = RT/P
    return vl,vv
end

function x0_sat_pure_hermite_spinodal(model,Vl00,Vv00,T)
    p(x) = pressure(model,x,T)
    TT = oneunit(eltype(model))*oneunit(T/T)
    Vl_spinodal,Vv_spinodal = Vl00*TT,Vv00*TT
    #=
    strategy:

    we create a quintic hermite interpolant for the pressure, given
    p,dpdv,d2pdv2 for the liquid and vapour phases.

    the hermite polynomial is translated in relation to Vl0: y = V - Vl0.
    we find yl,yv such as dpdv(yl) = dpdv(yv) = 0
    we recompute the new Vl0 and Vv0 and recalculate the hermite polynomial interpolant.
    we iterate until yl < tolerance.

    once the liquid spinodal is found, we can obtain the vapour one via optimization,
    using the current Vv estimate and the liquid spinodal as bounds.

    finally, once liquid and gas spinodals are found, we can obtain the spinodal initial point:
    psat ≈ 0.5*(p(liquid spinodal) + p(vapour spinodal))
    =#
    for i in 1:30
        fl,dfl,d2fl = Solvers.f∂f∂2f(p,Vl_spinodal)
        fv,dfv,d2fv = Solvers.f∂f∂2f(p,Vv_spinodal)
        poly_l = Solvers.hermite5_poly(Vl_spinodal,Vv_spinodal,fl,fv,dfl,dfv,d2fl,d2fv)
        dpoly_l = Solvers.polyder(poly_l)
        d2poly_l = Solvers.polyder(dpoly_l)
        function f0_l(x)
            df,d2f = evalpoly(x,dpoly_l),evalpoly(x,d2poly_l)
            return df,df/d2f
        end

        poly_v = Solvers.hermite5_poly(Vv_spinodal,Vl_spinodal,fv,fl,dfv,dfl,d2fv,d2fl)
        dpoly_v = Solvers.polyder(poly_v)
        d2poly_v = Solvers.polyder(dpoly_v)

        function f0_v(x)
            df,d2f = evalpoly(x,dpoly_v),evalpoly(x,d2poly_v)
            return df,df/d2f
        end

        prob_vl = Roots.ZeroProblem(f0_l,zero(Vl_spinodal))
        yl = Roots.solve(prob_vl,Roots.Newton())
        Vl_spinodal += yl
        #find proper bounds for yv
        ub_v = zero(Vv_spinodal)
        lb_v = Vl_spinodal - Vv_spinodal
        f_ub_v = evalpoly(zero(Vv_spinodal),dpoly_v)
        x = LinRange(lb_v,ub_v,10)

        #we need an initial point between the liquid spinodal and the initial point of the gas spinodal.
        for i in 1:9
            f_i = evalpoly(x[i],dpoly_v)
            if sign(f_i*f_ub_v) == -1
                lb_v = x[i]
                break
            end
            #cannot find any spinodal point. bail out.
            i == 9 && return Vl00,Vv00
        end
        yv_bounds = (lb_v,ub_v)
        prob_vv = Roots.ZeroProblem(f0_v,yv_bounds)
        yv = Roots.solve(prob_vv,Roots.LithBoonkkampIJzermanBracket())
        Vv_spinodal += yv

        if abs(yl) < sqrt(eps(Vl_spinodal)) && abs(yv) < sqrt(eps(Vv_spinodal))
            #yl found calculate mean spinodal pressure.
            P_max = evalpoly(yv,poly_v)
            P_min = evalpoly(yl,poly_l)
            P = 0.5*(P_max + P_min)
            Vv0 = volume(model,P,T,vol0 = Vv00)::typeof(TT)
            Vl0 = volume(model,P,T,vol0 = Vl00)::typeof(TT)
            #this initial point gives good estimated for the saturated vapour volume
            #but overpredicts the saturated liquid volume, causing failure in subsequent eq solvers.
            #isofugacity criteria does not work here.
            #On PCSAFT("water"), T = 0.9999Tc, this fails, even as the initial guesses provided
            #here are closer to the eq values than the ones sugested by x0_sat_pure_crit.
            return Vl0,Vv0
        end
    end
    return Vl00,Vv00 #failed to find spinodal
end

function vdw_x0_sat_pure(T,T_c,P_c,V_c)
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

#=
if we are near the critical point and we know it, we can use corresponding states to
obtain a good guess
=#

function _propaneref_rholsat(T)
    T_c = 369.89
    ρ_c = 5000.0
    T>T_c && return zero(T)/zero(T)
    Tr = T/T_c
    θ = 1.0-Tr
    ρ_l = (1.0 + 1.82205*θ^0.345 + 0.65802*θ^0.74 + 0.21109*θ^2.6 + 0.083973*θ^7.2)*ρ_c
    return ρ_l
end

function _propaneref_rhovsat(T)
    T_c = 369.89
    ρ_c = 5000.0
    T>T_c && return zero(T)/zero(T)
    Tr = T/T_c
    θ = 1.0 - Tr
    log_ρ_v_ρ_c = (-2.4887*θ^0.3785 -5.1069*θ^1.07 -12.174*θ^2.7 -30.495*θ^5.5 -52.192*θ^10 -134.89*θ^20)
    ρ_v = exp(log_ρ_v_ρ_c)*ρ_c
    return ρ_v
end

function x0_sat_pure_crit(model,T,crit::Tuple)
    Tc,Pc,Vc = crit
    return x0_sat_pure_crit(model,T,Tc,Pc,Vc)
end

function x0_sat_pure_crit(model,T,Tc,Pc,Vc)
    h = Vc*5000
    T0 = 369.89*T/Tc
    Vl0 = (1.0/_propaneref_rholsat(T0))*h
    Vv0 = (1.0/_propaneref_rhovsat(T0))*h
    _1 = SA[1.0]
    return Vl0,Vv0
end

function x0_sat_pure_crit(model,T)
    Tc,Pc,Vc = crit_pure(model)
    return x0_sat_pure_crit(model,T,Tc,Pc,Vc)
end

function scale_sat_pure(model)
    p    = 1/p_scale(model,SA[1.0])
    μ    = 1/Rgas(model)/T_scale(model,SA[1.0])
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
        px = exp(lnp̃)*ps
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

function x0_saturation_temperature(model,p)
    single_component_check(x0_saturation_temperature,model)
    coeffs = antoine_coef(model)
    if coeffs !== nothing
        return x0_saturation_temperature_antoine_coeff(model,p,coeffs)
    end
    #if obtaining the critical point is cheap, models can opt in by defining:
    #=
    x0_saturation_temperature(model::MyModel,p) = x0_saturation_temperature(model,p,crit_pure(model))
    =#
    return x0_saturation_temperature(model,p,nothing)
end

function x0_saturation_temperature(model::EoSModel,p,::Nothing)
    single_component_check(x0_saturation_temperature,model)
    return x0_saturation_temperature_refine(model,p)
end

function x0_saturation_temperature(model::EoSModel,p,crit::Tuple)
    single_component_check(x0_saturation_temperature,model)
    return x0_saturation_temperature_crit(model,p,crit)
end

#model has knowledge of the Antoine coefficients. use those to create an initial point.
function x0_saturation_temperature_antoine_coeff(model,p,coeffs)
    A,B,C = antoine_coef(model)
    lnp̄ = log(p / p_scale(model))
    T0 = T_scale(model)*(B/(A-lnp̄)-C)
    return x0_saturation_temperature_refine(model,p,T0)
end

function x0_saturation_temperature_crit(model::EoSModel,p,crit)
    Tc,Pc,Vc = crit
    A,B,C = (6.668322465137264,6.098791871032391,-0.08318016317721941) #universal antoine constants (RK)
    if Pc < p
        nan = zero(p)/zero(p)
        return (nan,nan,nan)
    end
    lnp̄ = log(p / Pc)
    T0 = Tc*(B/(A-lnp̄)-C)
    return x0_saturation_temperature_refine(model,p,T0)
end

#refine an initial temperature for x0_saturation_pressure, via calculating the saturation pressure at that temperature
#and performing second order extrapolation.
function dpdTsat_step(model,p,T0,crit::Union{Nothing,Tuple})
    #we want to use the crit point if available.
    return dpdTsat_step(model,p,T0,ChemPotVSaturation(;crit,crit_retry = isnothing(crit)))
end

function dpdTsat_step(model,p,T0)
    return dpdTsat_step(model,p,T0,ChemPotVSaturation(crit_retry = false))
end

function dpdTsat_step(model,p,T0,satmethod,multiple::Bool = true)
    T = T0*oneunit(eltype(model))*oneunit(p)*1.0
    dT = one(T)/zero(T)
    nan = zero(T)/zero(T)
    sat = (nan,nan,nan)
    n = multiple ? 10 : 1
    for i in 1:n
        sat = saturation_pressure(model,T,satmethod)
        pii,vli,vvi = sat
        if isnan(pii)
            return zero(pii)/zero(pii), sat
        end
        #=
        we use 1/T = 1/T0 + k*log(p/p0)
        k = -p0/(dpdT*T*T) 
        dpdT(saturation) = Δs/Δv (Clapeyron equation)
        =#
        dpdT = (VT_entropy(model,vvi,T) - VT_entropy(model,vli,T))/(vvi - vli)
        dTinvdlnp = -pii/(dpdT*T*T)
        Δlnp = log(p/pii)
        #dT = clamp(dTdp*Δp,-0.5*T,0.5*T)
        Tinv0 = 1/T
        Tinv = Tinv0 + dTinvdlnp*Δlnp
        dT = T - 1/Tinv
        T = 1/Tinv
        #!multiple && return T,sat
        if abs(dT)/T < 0.02
            return T,sat
        end
    end
    return T,sat
end

function x0_saturation_temperature_refine(model,p,T0::XX = 0.9*T_scale(model)*oneunit(p)*1.0,refine::Bool = true) where XX
    t,sat = dpdTsat_step(model,p,T0)
    _,vl,vv = sat
    vl = volume(model,p,t,phase = :liquid)
    vv = volume(model,p,t,phase = :gas)
    return t,vl,vv
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

"""
    solve_2ph_taylor(v10,v20,a1,da1,d2a1,a2,da2,d2a2,p_scale = 1.0,μ_scale = 1.0)

Solves the 2-phase problem with 1 component, using a 2nd order taylor aprox in helmholtz energy and a isothermal compressibility factor aproximation for pressure.
"""
function solve_2ph_taylor(v10,v20,a1,da1,d2a1,a2,da2,d2a2,p_scale = 1.0,μ_scale = 1.0)
    function F0(x)
        logv1,logv2 = x[1],x[2]
        v1,v2 = exp(logv1),exp(logv2)
        p1 = log(v1/v10)*(-v1*d2a1) - da1
        p2 = log(v2/v20)*(-v2*d2a2) - da2
        Δv1 = (v1 - v10)
        Δv2 = (v2 - v20)
        A1 = evalpoly(Δv1,(a1,da1,0.5*d2a1))
        A2 = evalpoly(Δv2,(a2,da2,0.5*d2a2))
        μ1 = A1 + p1*v1
        μ2 = A2 + p2*v2
        #F[1] = (μ1 - μ2)*μ_scale
        #F[2] = (p1 - p2)*p_scale
        F1 = (μ1 - μ2)*μ_scale
        F2 = (p1 - p2)*p_scale
        #return F
        return SVector((F1,F2))
    end
    x0 = SVector((log(v10),log(v20)))
    x = Solvers.nlsolve2(F0,x0,Solvers.Newton2Var())
    return exp(x[1]), exp(x[2])
end

function solve_2ph_taylor(model1::EoSModel,model2::EoSModel,T,v1,v2,p_scale = 1.0,μ_scale = 1.0)
    z = SA[1.0]
    f1(_V) = eos(model1,_V,T,z)
    f2(_V) = eos(model2,_V,T,z)
    a1,da1,d2a1 = Solvers.f∂f∂2f(f1,v1)
    a2,da2,d2a2 = Solvers.f∂f∂2f(f2,v2)
    return solve_2ph_taylor(v1,v2,a1,da1,d2a1,a2,da2,d2a2,p_scale,μ_scale)
end