"""
    x0_volume_liquid(model,T,z)
    x0_volume_liquid(model,p,T,z)

Returns an initial guess to the liquid volume, dependent on temperature `T` and composition `z`. By default is 1.25 times [`lb_volume`](@ref).
"""
function x0_volume_liquid(model,T,z)
    v_lb = lb_volume(model,T,z)
    return v_lb*1.25
end

x0_volume_liquid(model,p,T,z) = x0_volume_liquid(model,T,z)
x0_volume_liquid(model,T) = x0_volume_liquid(model,T,SA[1.0])
"""
    x0_volume_gas(model,p,T,z)

Returns an initial guess to the gas volume, depending of pressure `p`, temperature `T` and composition `z`. By default uses [`volume_virial`](@ref)
"""
function x0_volume_gas(model,p,T,z)
    B = second_virial_coefficient(model,T,z)
    nRT = sum(z)*Rgas(model)*T
    pmax = -0.25*nRT/B
    if B >= 0 || !isfinite(B)
        return nRT/p
    elseif pmax < p && B < 0
        return -2*B
    else
        return volume_virial(B,p,T,z)
    end
end

x0_volume_gas(model,p,T) = x0_volume_gas(model,p,T,SA[1.0])
"""
    x0_volume_solid(model,T,z)
    x0_volume_solid(model,p,T,z)

Returns an initial guess to the solid volume, dependent on temperature `T` and composition `z`. Needs to be defined for EoS that support solid phase. By default returns NaN. Can be overrided if the EoS defines `is_solid(::EoSModel) = true`
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
Returns `NaN` otherwise.
"""
function x0_volume(model, p, T, z = SA[1.0]; phase = :unknown)
    return x0_volume_impl(model,p,T,z,phase)
end

function x0_volume_impl(model, p, T, z = SA[1.0], phase = :unknown)
    if is_unknown(phase) || is_liquid(phase)
        return x0_volume_liquid(model,p,T,z)
    elseif is_vapour(phase)
        return x0_volume_gas(model,p,T,z)
    elseif is_supercritical(phase)
        return x0_volume_gas(model,p,T,z)
    elseif is_solid(phase)
        x0_volume_solid(model,p,T,z)
    else
        _0 = zero(Base.promote_eltype(model,p,T,z))
        return _0/_0
    end
end

"""
    lb_volume(model::EoSModel)
    lb_volume(model::EoSModel,z)
    lb_volume(model::EoSModel,T,z)

Returns the lower bound volume.
It has different meanings depending on the Equation of State, but symbolizes the minimum allowable volume at a certain composition:
- SAFT EoS: the packing volume
- Cubic EoS, covolume (b) parameter
On empiric equations of state, the value is chosen to match the volume of the conditions at maximum pressure and minimum temperature,
but the equation itself normally can be evaluated at lower volumes.
On SAFT and Cubic EoS, volumes lower than `lb_volume` will likely error.
The lower bound volume is used for guesses of liquid volumes at a certain pressure, saturated liquid volumes and critical volumes.

In most cases, the lower bound volume is independent of temperature. Some notable exceptions are the Quantum-Corrected Peng-Robinson cubic (`QCPR`) and Cubic-plus-Chain (CPC) models. For those,
it is better to define the three-argument variant `lb_volume(model,T,z)`
"""
function lb_volume end
lb_volume(model,T,z) = lb_volume(model,z)
lb_volume(model) = lb_volume(model,SA[1.0])
"""
    T_scale(model::EoSModel,z)
Represents a temperature scaling factor.

On any EoS based on Critical parameters (Cubic or Empiric EoS), the temperature scaling factor is chosen to be the critical temperature.
On SAFT or other molecular EoS, the temperature scaling factor is chosen to be a function of the potential depth ϵ.
Used as scaling factors in [`saturation_pressure`](@ref) and as input for solving [`crit_pure`](@ref)
"""
function T_scale end
T_scale(model) = T_scale(model,SA[1.0])
"""
    p_scale(model::EoSModel,z)
Represents a pressure scaling factor.

On any EoS based on Critical parameters (Cubic or Empiric EoS), the pressure scaling factor is
chosen to be a function of the critical pressure. On SAFT or other molecular EoS, the pressure scaling factor is chosen to a function of ∑(zᵢ*ϵᵢ*(σᵢᵢ)³)
Used as scaling factors in [`saturation_pressure`](@ref) and as input for solving [`crit_pure`](@ref).

By default, it can be defined as a function of `Clapeyron.lb_volume` and `Clapeyron.T_scale`
"""
function p_scale end

p_scale(model) = p_scale(model,SA[1.0])

function p_scale(model,z)
    Ts = T_scale(model,z)
    sum(z)*Rgas(model)*Ts/lb_volume(model,Ts,z)
end
"""
    antoine_coef(model)
Should return a 3-Tuple containing reduced Antoine Coefficients. The Coefficients follow the correlation:
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
    saturation_model(model) = model

Returns the model used to calculate pure saturation properties and pure critical points. By default returns the input model.

"""
saturation_model(model::T) where T = model

"""
    has_fast_crit_pure(model::EoSModel)::Bool

Used to indicate if a model can calculate their critical point without iterative calculations.
Having a critical point available results in speed ups for saturation calculations.
By default returns `false`.
"""
function has_fast_crit_pure(model)::Bool
    satmodel = saturation_model(model)
    if satmodel !== model
        return has_fast_crit_pure(satmodel)
    else
        return false
    end
end

"""
    x0_sat_pure(model::EoSModel,T)
    x0_sat_pure(model,T,crit)

Returns a 2-tuple corresponding to `(Vₗ,Vᵥ)`, where `Vₗ` and `Vᵥ` are the liquid and vapor initial guesses.
Used in [`saturation_pressure`](@ref) methods that require initial volume guesses.
It can be overloaded to provide more accurate estimates if necessary. If an EoS model provides a fast method for `crit_pure`, overloading `has_fast_crit_pure` will provide `x0_sat_pure` with additional information to improve its accuracy.
"""
function x0_sat_pure(model,T)
    satmodel = saturation_model(model)
    if satmodel !== model
        return x0_sat_pure(satmodel,T)
    end
    if !has_fast_crit_pure(model)
    _,vl,vv = x0_sat_pure_virial(model,T)
    else
    _,vl,vv = x0_sat_pure_crit(model,T)
    end
    return vl,vv
end

function x0_sat_pure(model,T,crit)
    satmodel = saturation_model(model)
    if satmodel !== model
        return x0_sat_pure(satmodel,T,crit)
    end

    if has_fast_crit_pure(model)
        #if the model has a custom method here, it will be dispatched to that one.
        return x0_sat_pure(model,T)
    end

    if isnothing(crit)
        _,vl,vv = x0_sat_pure_virial(model,T)
    else
        _,vl,vv = x0_sat_pure_crit(model,T,crit)
    end
    return vl,vv
end


"""
    p,vl,vv = x0_sat_pure_virial(model,T)

Calculates initial points for pure saturation pressure using a virial + corresponding states approach.
The corresponding states model (a vdW fluid fitted from 2 p-V points) is used to select between zero-pressure or spinodal initial points.
The points selected to fit the vdW fluid are a function of B(T).
"""
function x0_sat_pure_virial(model,T)
    z = SA[1.0]
    single_component_check(x0_sat_pure,model)

    #=theory as follows
    #given T = Teos:
    - calculate B(Teos),Vv = -2B,Pv = P(V = -2B)
    - Pl = f(B), Vl = volume(p = Pl)
    - given Pl,Vl,Pv,Vv, calculate vdW parameters a,b
    with a,b, calculate T̃ = RT*b/a. we use this parameter to determine if we are near the critical point or not.
    and switch to the appropiate strategy.
    =#

    R̄ = Rgas(model)
    RT= R̄*T
    B = second_virial_coefficient(model,T,z)
    _0,_1 = zero(B),oneunit(B)
    lb_v = lb_volume(model,T,z)*_1

    #=
    some very complicated models, like DAPT, fail on the calculation of the second virial coefficient.
    while this is a numerical problem in the model itself, it is better to catch this early.

    while a volume calculation is harder
    =#
    if isnan(B)
        x0l = 3*lb_v
        px = pressure(model,x0l,T)
        if px < 0 #low pressure
            p = RT/vv
            return p,vl,vv
        else #high pressure?
            vl = 4*lb_v
            return pressure(model,vl,T),vl,20*lb_v
        end
    end


    _0 = zero(B)
    #virial volume below lower bound volume.
    #that means that we are way over the critical point
    if -2B < lb_v
        _nan = _0/_0
        return (_nan,_nan,_nan)
    end
    vv_virial = -2*B
    pv_eos = pressure(model,vv_virial,T) #eos predicted pressure, gas phase
    #calculate a suitable liquid pressure from the virial coefficient.
    pl0 = liquid_pressure_from_virial(model,T,B,pv_eos)
    vl = volume(model,pl0,T,z,phase = :l)
    #=the basis is that p = RT/v-b - a/v2
    we can interpolate a vdW EoS between a liquid and a gas (p,v) point.
    with that, we solve for a and b
    with a and b, we calculate T̃ = RTb/a = Tr/α(Tr)
    at Tr = 1, T̃max is Ωb/Ωa. we can check T̃/T̃max to see how close (or far) we are to the critical point
    =#

    a,b = vdw_coefficients(vl,pl0,vv_virial,pv_eos,T)
    if b < lb_v
        b = lb_v
        a = RT*(b - B)
    end
    T̃ = RT*b/a
    T̃max = 0.2962962962962963 # Ωb/Ωa = 8/27
    T̃min = 0.1*T̃max
    if isnan(T̃)
        #fails on two ocassions:
        #near critical point, or too low.
        #return high pressure estimate
        return pressure(model,4*lb_v,T),4*lb_v,vv_virial + 2*lb_v
    elseif T̃ >= T̃max
        x0l = 4*lb_v
        x0v = vv_virial + 2*lb_v
        return pressure(model,4*lb_v,T),x0l,x0v
    elseif T̃ < T̃min
        #gas volume over threshold. but not diverged.
        #normally this happens at low temperatures. we could suppose that Vl0 is a
        #"zero-pressure" volume, apply corresponding strategy
        return x0_sat_pure_near0(model,T,vl;B = B)
    else
        psat,_,vv_B = x0_sat_pure_near0(model,T,vl,B = B,refine_vl = false)
        if T̃/T̃max > 0.55
            B_vdw = b - a/RT
            vv_vdw_2b = -2*B_vdw
            if vv_vdw_2b < vv_B
                ps,vls,vvs = x0_sat_pure_spinodal(model,T,vl,vv_B,B)
                return ps,vls,vvs
            end
        end
        return psat,vl,vv_B
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

"""
    p,vl,vv = x0_sat_pure_lk(model,T,crit,ω)

Calculates initial points for pure saturation pressure, using the Lee-Kesler correlation.
"""
function x0_sat_pure_lk(model,T,crit,ω)
    _0 = zero(Base.promote_eltype(model,T,ω))
    nan = _0/_0
    tc,pc,vc = crit
    T > tc && (return nan,nan,nan)
    tr = T/tc
    trinv = inv(tr)
    lntr = log(tr)
    tr6 = tr^6
    f0 = 5.92714 - 6.09648*trinv - 1.28862*lntr + 0.169347*tr6
    f1 = 15.2518 - 15.6875*trinv - 13.4721*lntr + 0.43577*tr6
    lnpr = f0 + ω*f1
    psat = exp(lnpr)*pc
    vl = volume(model,psat,T,phase = :l)
    vv = volume(model,psat,T,phase = :v)
    return psat,vl,vv
end

"""
    p,vl,vv = x0_sat_pure_near0(model,T,vl0 = volume(model,zero(T),T,phase = :l);
                                B = second_virial_coefficient(model,T),
                                refine_vl = true)

Calculates initial points for pure saturation pressure, using a zero-pressure volume approach.
If `refine_vl` is set to `true`, then the liquid volume will be recalculated using the calculated saturation pressure, otherwise it will be returned as is.
"""
function x0_sat_pure_near0(model, T, vl0 = volume(model,zero(T),T,phase = :l);B = second_virial_coefficient(model,T), refine_vl = true)
    R̄ = Rgas(model)
    z = SA[1.0]
    RT = R̄*T
    ares = a_res(model, vl0, T, z)
    lnϕ_liq0 = ares - 1 + log(RT/vl0)
    p = exp(lnϕ_liq0)
    pB = -0.25*RT/B
    vv = volume_virial(B,p,T)
    if refine_vl && pB/p > 10
        vl = volume(model,p,T,z,vol0 = vl0,phase = :l)
        if vl ≈ vv #refinement failed, stick with vl0
            vl = vl0*oneunit(vv)
        end
    else
        vl = vl0*oneunit(vv)
    end
    if isnan(vv)
        vv = RT/p
    end
    return p,vl,vv
end

function liquid_pressure_from_virial(model,T,B = second_virial_coefficient(model,T),pv_eos = pressure(model,-2*B,T))
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

    vv_virial = -2*B #maximum gas volume predicted by virial equation
    pv_virial = -0.25*Rgas(model)*T/B #maximum virial predicted pressure
    γT = pv_eos/pv_virial

    #this handles pv_eos = NaN and pv_eos < pv_virial, returning an equivalent result to using pv_eos = pv_virial
    !(pv_eos > pv_virial) && (return 1.12491990759086*pv_virial*oneunit(pv_eos))
    #fitted function, using all coolprop fluids, at Tr = 1
    aγ,bγ,cγ = 1.2442071971165476e-5, -8.695786307570637, 1.0505452946870144
    γc = aγ*exp(-γT*bγ) + cγ
    return γc*pv_virial #pressure of which we are (almost) sure there exists a liquid root
end

function pure_spinodal(model,_T,z = SA[1.0];phase = :l)
    B = second_virial_coefficient(model,_T)
    _v_ub = -2*B
    pv_eos = pressure(model,_v_ub,_T,z)
    pl = liquid_pressure_from_virial(model,_T,B,pv_eos)
    _v_lb = volume(model,pl,_T,z,phase = :l)
    T,v_lb,v_ub = promote(_T,_v_lb,_v_ub)
    return pure_spinodal(model,T,v_lb,v_ub,phase,true,z)
end

#given an hermite polynomial that interpolates the spinodals
#find the intermediate point where d2pdv2 = 0 and d3pdv3 < 0
#this point is in between the liquid and vapour spinodals.
function _find_vm(dpoly,v_lb::K,v_ub::K) where K
    d2poly = Solvers.polyder(dpoly)
    d3poly = Solvers.polyder(d2poly)
    lb = zero(v_lb)
    ub = v_ub - v_lb
    if iszero(last(d2poly))
        c,b,a,_ = d2poly
        if iszero(a)
            #bx + c = 0
            v = -c/b
            nr,v1,v2,v3 = 1,v,v,v
        else
            dd = sqrt(b*b - 4*a*c)
            isnan(dd) && return zero(K)/zero(K)
            v1 = (-b + dd)/(2*a)
            v2 = (-b - dd)/(2*a)
            v3 = zero(K)/zero(K)
            nr = 2
        end
    else
        nr,v1,v2,v3 = Solvers.real_roots3(d2poly)
    end

    if evalpoly(v1,dpoly) > 0 && (lb <= v1 <= ub) && evalpoly(v1,d3poly) < 0
        return v1 + v_lb
    elseif evalpoly(v2,dpoly) > 0 && (lb <= v2 <= ub) && nr > 1 && evalpoly(v2,d3poly) < 0
        return v2 + v_lb
    elseif evalpoly(v2,dpoly) > 0 && (lb <= v3 <= ub) && nr > 2 && evalpoly(v3,d3poly) < 0
        return v3 + v_lb
    else
        return zero(K)/zero(K)
    end
end

function pure_spinodal_newton_bracket(model,T,v,f,dp_scale,z = SA[1.0])
    vlo,vhi = v
    flo,fhi = f
    vs = 0.5*(vlo + vhi)
    atol = 1e-8
    vs_old = vs*Inf
    for j in 1:25
        pj,dpj,d2pj = p∂p∂2p(model,vs,T,z)
        Δ = dpj/d2pj
        vs_newton = vs - Δ
        fs = dpj
        if fs*flo < 0
            vhi = vs
            fhi = fs
        elseif fs*fhi < 0
            vlo = vs
            flo = fs
        else
            return vs
        end

        if vlo <= vs_newton <= vhi || abs(Δ)/vs < 0.01
            vs_old = vs
            vs = vs_newton
        else
            Δ = vs_old - 0.5*(vlo + vhi)
            vs_old = vs
            vs = 0.5*(vlo + vhi)

        end
        if abs(Δ) < atol || abs(dp_scale*fs) < atol
            return vs
        end
    end

    return zero(vs)/zero(vs)
end

function pure_spinodal_newton(model,T,z,v0,dp_scale = v0*v0/(Rgas(model)*T))
    function dp(vs) #dpdrho = 0
        p(rho) = pressure(model,1/rho,T,z)
        pj,dpj,d2pj = Solvers.f∂f∂2f(p,1/vs)
        return dpj/dp_scale,dpj/d2pj
    end

    prob = Roots.ZeroProblem(dp,1/v0)
    v = Roots.solve(prob,Roots.Newton())
end

function pure_spinodal(model,T::K,v_lb::K,v_ub::K,phase::Symbol,retry,z = SA[1.0]) where K
    fl,dfl,d2fl = p∂p∂2p(model,v_lb,T,z)
    fv,dfv,d2fv = p∂p∂2p(model,v_ub,T,z)
    dfx = ifelse(is_liquid(phase),dfl,dfv)
    vx = ifelse(is_liquid(phase),v_lb,v_ub)
    nan = zero(fl)/zero(fl)
    isnan(vx) && return nan
    poly = Solvers.hermite5_poly(v_lb,v_ub,fl,fv,dfl,dfv,d2fl,d2fv)
    dpoly = Solvers.polyder(poly)

    dp_scale = evalpoly(vx - v_lb,dpoly)
    #we already have a bracket.
    if dfl*dfv < 0
        return pure_spinodal_newton_bracket(model,T,(v_lb,v_ub),(dfl,dfv),dp_scale,z)
    end

    #find the middle point between the liquid and vapour spinodals.
    vm = _find_vm(dpoly,v_lb,v_ub)
    fm,dfm,d2fm = p∂p∂2p(model,vm,T,z)
    #find the liquid of gas spinodal using the quintic hermite interpolation.
    v_bracket_hermite = minmax(vx - v_lb,vm - v_lb)
    !(evalpoly(vx - v_lb,dpoly)*evalpoly(vm - v_lb,dpoly) < 0) && return nan
    v_spinodal_hermite_prob = Roots.ZeroProblem(Base.Fix2(evalpoly,dpoly),v_bracket_hermite)
    vh = Roots.solve(v_spinodal_hermite_prob,xrtol = 1e-5) + v_lb
    fh,dfh,d2fh = p∂p∂2p(model,vh,T,z)
    unstable_not_found = dfx < 0 && dfm < 0 && dfh < 0

    if unstable_not_found
        if !retry
            return nan
        end

        v_ub_new = v_ub
        v_lb_new = v_lb - dfl/d2fl
        if dfm < 0 && d2fm < 0 && vm < v_ub_new
            v_ub_new = vm
        end

        if dfh < 0 && d2fh < 0 && vh < v_ub_new
            v_ub_new = vh
        end

        phase_h = VT_identify_phase(model,vh,T,z)

        if is_vapour(phase_h) && is_liquid(phase) && d2fh > 0 && d2fm > 0
            #v_lb_new = v_lb - dfl/d2fl
            v_ub_new = vh
        end

        if is_liquid(phase_h) && is_liquid(phase) && vh > v_lb_new
            v_lb_new = vh
        end

        return pure_spinodal(model,T,v_lb_new,v_ub_new,phase,false,z)
    end

    if dfx*dfh <= 0
        if vx < vh
            v_bracket = (vx,vh)
            dp_bracket = (dfx,dfh)
        else
            v_bracket = (vh,vx)
            dp_bracket = (dfh,dfx)
        end
        return pure_spinodal_newton_bracket(model,T,v_bracket,dp_bracket,dp_scale,z)
    elseif dfx*dfm <= 0
        if vx < vm
            v_bracket = (vx,vm)
            dp_bracket = (dfx,dfm)
        else
            v_bracket = (vm,vx)
            dp_bracket = (dfm,dfx)
        end
        return pure_spinodal_newton_bracket(model,T,v_bracket,dp_bracket,dp_scale,z)
    else
        throw(error("Cannot determine spinodal bracket for $(typeof(model)) at phase = :$phase. input volume values are: ($v_lb,$v_ub)"))
    end
end

"""
    p,vl,vv = x0_sat_pure_spinodal(model,T,B = second_virial_coefficient(model,T))

Calculates initial points for pure saturation pressure, using a spinodal approach.
The saturation pressure is assumed to be `(psl + psv)/2` where `psl` and `psv` are the pressures of the liquid and vapour spinodals.
"""
function x0_sat_pure_spinodal(model,T,B = second_virial_coefficient(model,T))
    v_ub = -2*second_virial_coefficient(model,T)
    pl = liquid_pressure_from_virial(model,T,B)
    v_lb = volume(model,pl,T,phase = :l)
    return x0_sat_pure_spinodal(model,T,v_lb,v_ub,B)
end

function x0_sat_pure_spinodal(model,T,v_lb,v_ub,B = second_virial_coefficient(model,T),Vc = nothing)
    if Vc === nothing
        vc = zero(v_lb)/zero(v_ub)
    else
        vc,_,_ = promote(Vc,v_lb,v_ub)
    end

    p(x) = pressure(model,x,T)

    if isnan(vc)
        vsl = pure_spinodal(model,T,v_lb,v_ub,:l,true)
    else
        vsl = pure_spinodal(model,T,v_lb,vc,:l,true)
    end

    psl = p(vsl)
    if isnan(vsl)
        return pressure(model,v_lb,T),v_lb,v_ub
    end
    if isnan(vc)
        vsv = pure_spinodal(model,T,vsl,v_ub,:v,true)
    else
        vsv = pure_spinodal(model,T,vc,v_ub,:v,true)
    end

    if isnan(vsv)
        return pressure(model,v_lb,T),v_lb,v_ub
    end

    plb = p(v_lb)
    pub = p(v_ub)
    psl = p(vsl)
    psv = p(vsv)
    pmid = 0.5*max(zero(psl),psl) + 0.5*psv
    #=
    vl = volume(model,pmid,T,phase = :l,vol0 = v_lb)
    vv = volume(model,pmid,T,phase = :v)
    return pmid,vl,vv
    =#
    #
    if plb < psv
        vsl_lb = volume(model,psv,T,phase = :l, vol0 = v_lb)
    else
        vsl_lb = one(psl)*v_lb
    end

    if pub >= min(pmid,max(psl,zero(psl)))
        vsv_ub = Rgas(model)*T/pmid
    else
        vsv_ub = one(psv)*v_ub
    end
    return _x0_sat_pure_spinodal(model,T,vsl_lb,vsv_ub,vsl,vsv,B)
end

function _x0_sat_pure_spinodal(model,T,vsl_lb,vsv_ub,vsl,vsv,B)
    psl,_,d2psl = p∂p∂2p(model,vsl,T,SA[1.0])
    psv,_,d2psv = p∂p∂2p(model,vsv,T,SA[1.0])
    psl_lb,dpsl_lb,d2psl_lb = p∂p∂2p(model,vsl_lb,T,SA[1.0])
    dpsl = zero(psl)
    poly_l = Solvers.hermite5_poly(vsl_lb,vsl,psl_lb,psl,dpsl_lb,dpsl,d2psl_lb,d2psl)
    ps_mid = 0.5*(psv + max(psl,zero(psl)))
    vl = volume_from_spinodal(ps_mid,poly_l,vsl_lb,0.5*(vsl_lb + vsl) - vsl_lb)
    if psl < 0
        vv = volume_virial(B,ps_mid,T)
        isnan(vv) && (vv = Rgas(model)*T/ps_mid)
        return ps_mid,vl,vv
    end
    psv_ub,dpsv_ub,d2psv_ub = p∂p∂2p(model,vsv_ub,T,SA[1.0])
    dpsv = zero(psl)
    poly_v = Solvers.hermite5_poly(vsv,vsv_ub,psv,psv_ub,dpsv,dpsv_ub,d2psv,d2psv_ub)
    vv = volume_from_spinodal(ps_mid,poly_v,vsv,(zero(vsv),vsv_ub - vsv))
    return ps_mid,vl,vv
end

function volume_from_spinodal(p,poly,vshift,v0)
    f(v) = p - evalpoly(v,poly)


    if length(v0) == 2
        v1,v2 = v0
        f1,f2 = f(v1),f(v2)
        if f1*f2 < 0
            prob = Roots.ZeroProblem(f,v0)
            return Roots.solve(prob) + vshift
        else
            #something really wrong happened with the bracketing,
            #hopefully the hermite polynomial should reproduce the EoS in a vicinity of the
            #interpolated region.
            if abs(f1) < abs(f2)
                prob = Roots.ZeroProblem(f,v1)
                return Roots.solve(prob) + vshift
            else
                prob = Roots.ZeroProblem(f,v2)
                return Roots.solve(prob) + vshift
            end
        end
    end
    prob = Roots.ZeroProblem(f,v0)
    return Roots.solve(prob) + vshift
end

x0_sat_pure_crit(model,T) = x0_sat_pure_crit(model,T,crit_pure(model))

"""
    p,vl,vv = x0_sat_pure_crit(model,T,crit = crit_pure(model))

Calculates initial points for pure saturation pressure, using a combination of spinodal + zero-pressure + critical extrapolation approaches.
The strategy selected depends on the value of `T/Tc`
"""
function x0_sat_pure_crit(model,_T,crit::NTuple{3,Any})
    _Tc,_Pc,_Vc = crit
    _1 = one(Base.promote_eltype(model,_T,_Tc,_Pc,_Vc))
    _0 = zero(_1)
    _,T,Tc,Pc,Vc = promote(_1,_T,_Tc,_Pc,_Vc)
    Tr = T/Tc
    nan = _0/_0
    RT = Rgas(model)*T
    if Tr == 1
        return Pc,Vc,Vc
    elseif Tr > 1
        return nan,nan,nan
    elseif 0.99 < Tr < 1.0
        vl,vv = critical_vsat_extrapolation(model,T,Tc,Vc)
        p = pressure(model,vl,T)
        return p,vl,vv
    end

    B = second_virial_coefficient(model,T)
    v_ub = -2B

    if v_ub < 0
        return x0_sat_pure_near0(model,T,B = B)
    end

    pl0 = liquid_pressure_from_virial(model,T,B)
    v_lb = volume(model,pl0,T,phase = :l)

    if 0.8 <= Tr <= 0.99
        return x0_sat_pure_spinodal(model,T,v_lb,v_ub,B,Vc)
    elseif 0 <= Tr < 0.8
        return x0_sat_pure_near0(model,T,v_lb,B = B)
    else
        return nan,nan,nan
    end
end

function equilibria_scale(model,z = SA[1.0])
    Ts = T_scale(model,z)
    lb_v = lb_volume(model,Ts,z)
    R = Rgas(model)
    μ = 1/(R*Ts)
    p = lb_v*μ
    #p    = 1/p_scale(model,SA[1.0])
    #μ    = 1/Rgas(model)/T_scale(model,SA[1.0])
    return p,μ
end

"""
    x0_psat(model::EoSModel, T,crit = nothing)
Initial point for saturation pressure, given the temperature and V,T critical coordinates.
On moderate pressures it will use a Zero Pressure initialization. On pressures near the critical point it will switch to spinodal finding.
Used in [`saturation_pressure`](@ref) methods that require initial pressure guesses.
If the initial temperature is over the critical point, it returns `NaN`.
It can be overloaded to provide more accurate estimates if necessary.
"""
function x0_psat(model,T)
    satmodel = saturation_model(model)
    satmodel !== model && return x0_psat(satmodel,T)
    if has_fast_crit_pure(model)
        p,_,_ = x0_sat_pure_virial(model,T)
    else
        p,_,_ = x0_sat_pure_crit(model,T)
    end
    return p
end

function x0_psat(model,T,crit)
    satmodel = saturation_model(model)
    satmodel !== model && return x0_psat(satmodel,T,crit)
    has_fast_crit_pure(model) && return x0_psat(model,T)
    if isnothing(crit)
        return first(x0_sat_pure_virial(model,T))
    else
        return first(x0_sat_pure_crit(model,T,crit))
    end
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
    #@show coeffs
    #if coeffs !== nothing
    #    return x0_saturation_temperature_antoine_coeff(model,p,coeffs)
    #end
    #if obtaining the critical point is cheap, models can opt in by defining:
    #=
    x0_saturation_temperature(model::MyModel,p) = x0_saturation_temperature(model,p,crit_pure(model))
    =#
    if !has_fast_crit_pure(model)
        return x0_saturation_temperature_refine(model,p)
    else
        return x0_saturation_temperature_crit(model,p,crit_pure(model))
    end
end

function x0_saturation_temperature(model::EoSModel,p,::Nothing)
    single_component_check(x0_saturation_temperature,model)
    if has_fast_crit_pure(model)
        return x0_saturation_temperature_crit(model,p,crit_pure(model))
    else
        return x0_saturation_temperature_refine(model,p)
    end
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
    T0 = critical_tsat_extrapolation(model,p,crit)
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
        dpdT = dpdT_saturation(model,vli,vvi,T)
        dTinvdlnp = -pii/(dpdT*T*T)
        Δlnp = log(p/pii)
        #dT = clamp(dTdp*Δp,-0.5*T,0.5*T)
        Tinv0 = 1/T
        Tinv = Tinv0 + dTinvdlnp*Δlnp
        dT = T - 1/Tinv
        if 1/Tinv > T
            T = 0.5*T + 0.5/Tinv #we could skip over the critical temperature
        else
            T = 1/Tinv
        end
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
    x0_crit_pure(model::EoSModel,z)
Returns a 2-tuple corresponding to
    `(k,log10(Vc0))`, where `k` is `Tc0/T_scale(model,z)`
"""
function x0_crit_pure end

x0_crit_pure(model) = x0_crit_pure(model,SA[1.0])

function x0_crit_pure(model::EoSModel,z)
    Ts = T_scale(model,z)
    lb_v = lb_volume(model,Ts,z)/sum(z)
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

Solves the 2-phase problem with 1 component, using a 2nd order taylor approx in Helmholtz energy and a isothermal compressibility factor approximation for pressure.
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
    v1,v2 = exp(x[1]), exp(x[2])
    p1 = log(v1/v10)*(-v1*d2a1) - da1
    return v1, v2, p1
end

function solve_2ph_taylor(model1::EoSModel,model2::EoSModel,T,v1,v2,p_scale = 1.0,μ_scale = 1.0)
    z = SA[1.0]
    f1(_V) = eos(model1,_V,T,z)
    f2(_V) = eos(model2,_V,T,z)
    a1,da1,d2a1 = Solvers.f∂f∂2f(f1,v1)
    a2,da2,d2a2 = Solvers.f∂f∂2f(f2,v2)
    return solve_2ph_taylor(v1,v2,a1,da1,d2a1,a2,da2,d2a2,p_scale,μ_scale)
end

"""
    critical_vsat_extrapolation(model,T,Tc,Vc)
    critical_vsat_extrapolation(model,T,crit = crit_pure(model))

Given critical information and a temperature, extrapolate the liquid and vapor saturation volumes.
"""
function critical_vsat_extrapolation(model,T,Tc,Vc,z = SA[1.0])
    if T > Tc
        _0 = zero(Base.promote_eltype(model,T,z))
        nan = _0/_0
        return nan,nan
    end
    ρc = 1/Vc
    function dp(ρ,T)
        _,dpdV = p∂p∂V(model,1/ρ,T,z)
        return -sum(z)*dpdV*ρ*ρ
    end
    #Solvers.derivative(dρ -> pressure(model, 1/dρ, T), ρ)
    _,d2p,d3p = Solvers.∂J2(dp,ρc,Tc)
    ∂²p∂ρ∂T = d2p[2]
    ∂³p∂ρ³ = d3p[1,1]
    Bp = sqrt(6 * Tc * ∂²p∂ρ∂T / ∂³p∂ρ³)
    ΔT = (Tc - T)/T
    Δρ = Bp*sqrt(ΔT)
    ρl = ρc + Δρ
    ρv = ρc - Δρ
    return 1/ρl,1/ρv
end

critical_vsat_extrapolation(model,T) = critical_vsat_extrapolation(model,T,crit_pure(model))
critical_vsat_extrapolation(model,T,crit) = critical_vsat_extrapolation(model,T,crit[1],crit[3])

"""
    critical_psat_extrapolation(model,T,Tc,Pc,Vc)
    critical_psat_extrapolation(model,T,Tc,Vc)
    critical_psat_extrapolation(model,T,crit = crit_pure(model))

Given critical information and a temperature, extrapolate the saturation pressure.

!!! note
    This function will not check if the input temperature is over the critical point.
"""
function critical_psat_extrapolation(model,T,Tc,Pc,Vc)
    _p(_T) = pressure(model,Vc,_T)
    dpdT = Solvers.derivative(_p,Tc)
    dTinvdlnp = -Pc/(dpdT*Tc*Tc)
    Δlnp = (1/T - 1/Tc)/dTinvdlnp
    p = exp(Δlnp)*Pc
end

critical_psat_extrapolation(model,T) = critical_psat_extrapolation(model,T,crit_pure(model))
critical_psat_extrapolation(model,T,crit) = critical_psat_extrapolation(model,T,crit[1],crit[2],crit[3])
critical_psat_extrapolation(model,T,Tc,Vc) = critical_psat_extrapolation(model,T,Tc,pressure(model,Vc,Tc),Vc)


"""
    critical_tsat_extrapolation(model,p,Tc,Pc,Vc)
    critical_tsat_extrapolation(model,p,Tc,Vc)
    critical_tsat_extrapolation(model,p,crit = crit_pure(model))

Given critical information and a pressure, extrapolate the saturation temperature.

!!! note
    This function will not check if the input pressure is over the critical point.

"""
function critical_tsat_extrapolation(model,p,Tc,Pc,Vc,z = SA[1.0])
    _p(_T) = pressure(model,Vc,_T,z)
    dpdT = Solvers.derivative(_p,Tc)
    dTinvdlnp = -Pc/(dpdT*Tc*Tc)
    Δlnp = log(p/Pc)
    Tinv = 1/Tc + dTinvdlnp*Δlnp
    T = 1/Tinv
end

critical_tsat_extrapolation(model,p) = critical_tsat_extrapolation(model,p,crit_pure(model))
critical_tsat_extrapolation(model,p,crit) = critical_tsat_extrapolation(model,p,crit[1],crit[2],crit[3])
critical_tsat_extrapolation(model,p,Tc,Vc) = critical_tsat_extrapolation(model,p,Tc,pressure(model,Vc,Tc),Vc)


dpdT_saturation(model::EoSModel,v1::Number,v2,T) = dpdT_saturation(model,model,v1,v2,T,SA[1.0],SA[1.0])
dpdT_saturation(model1::EoSModel,model2::EoSModel,v1,v2,T) = dpdT_saturation(model1,model2,v1,v2,T,SA[1.0],SA[1.0])
function dpdT_saturation(model1::EoSModel,model2::EoSModel,v1,v2,T,w1,w2)
    ∑w1 = sum(w1)
    ∑w2 = sum(w2)

    dS_res = VT_entropy_res(model1,v1,T,w1)/∑w1 - VT_entropy_res(model2,v2,T,w2)/∑w2

    R1,R2 = Rgas(model1),Rgas(model2)
    ∑fx1,∑fx2 = R1*sum(xlogx,w1)/∑w1,R2*sum(xlogx,w2)/∑w2
    Δx = ∑fx1 - ∑fx2
    ∑fv1,∑fv2 = R1*∑w1*log(∑w1*v1), R2*∑w2*log(∑w1*v2)
    Δv = ∑fv1 - ∑fv2
    dS_ideal = Δx + Δv #Rgas(model1)*(log(v1/v2)
    dS = dS_res + dS_ideal
    dv = (v1/∑w1 - v2/∑w2)
    return dS/dv
end
