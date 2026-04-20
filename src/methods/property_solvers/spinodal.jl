
"""
    spinodal_pressure(model::EoSModel, T, x; v0, phase)

Calculates the spinodal pressure `p` and `V` volume for a given temperature `T` and composition `z`. Returns a tuple, containing:
- Spinodal pressure `[Pa]`
- Spinodal volume `[m³]`

Calculates either the liquid or the vapor spinodal point depending on the given starting volume `v0` or the `phase`. The keyword `phase` is ignored if `v0` is given.
"""
function spinodal_pressure(model::EoSModel,T,z=SA[1.];v0=nothing,phase=:unknown)
    check_arraysize(model,z)

    is_liquid(phase) || is_vapour(phase) || !isnothing(v0) || error("Either `v0` or `phase` has to be specified!")

    ∑z = sum(z)
    if isnothing(v0) && !is_unknown(phase) && length(model) == 1
        V_spin = pure_spinodal(model,T;phase)
        return pressure(model,V_spin,T),V_spin*∑z
    end

    model, idx_r = index_reduction(model,z)
    x = z[idx_r]

    #reduced model has length = 1
    if isnothing(v0) && !is_unknown(phase) && length(model) == 1
        V_spin = pure_spinodal(model,T;phase)
        return pressure(model,V_spin,T),V_spin*∑z
    end

    # Determine initial guess (if not provided)
    if !isnothing(v0)
        _v0 = one(Base.promote_eltype(model,T,z))*v0
    else
        _v0 = x0_spinodal_pressure(model,T,z,phase)
    end

    f(vz) = det_∂²A∂ϱᵢ²(model, T, x ./ exp(vz))
    if !isfinite(_v0)
        nan = zero(_v0)/zero(_v0)
        return nan,nan
    end
    prob = Roots.ZeroProblem(Solvers.to_newton(f),log(_v0))
    log_V_spin = Roots.solve(prob,Roots.Newton())
    V_spin = exp(log_V_spin)
    p_spin = pressure(model,V_spin,T,x)
    return p_spin,V_spin
end

function x0_spinodal_pressure(model, T, z = SA[1.0], phase = :unknown)

    is_liquid(phase) || is_vapour(phase) || error("phase` has to be specified!")

    v0_edge,pure_sats = x0_edge_pressure(model,T,z)
    p0_bubble,p0_dew = v0_edge
    pmin_sat,pmax_sat = extrema(first,pure_sats)
    edge,crit,status = _edge_pressure(model,T,z,v0_edge)
    P_edge,v_l,v_v = edge

    TT = Base.promote_eltype(model,P_edge,T,z)
    xx = similar(z,TT)

    if status == :failure
        return zero(TT)/zero(TT)
    end

    vx = if status == :supercritical
        crit[3]
    else
        if is_liquid(phase)
            v_l
        else
            v_v
        end
    end

    xx .= z ./ vx

    detx = det_∂²A∂ϱᵢ²(model,T,xx)

    status != :supercritical && detx > 0 && (return TT(vx))

    if status == :supercritical && detx > 0
        vmin = volume(model,pmax_sat,T,z,phase = :l)
        vmax = volume_virial(model,pmin_sat,T,z)
        f(vvz) = let model = model,T = T,z = z,xx  = xx
            xx .= z ./ vvz
            det_∂²A∂ϱᵢ²(model,T,xx)::TT
        end
        res = Solvers.optimize(f,(vmin,vmax),Solvers.BoundOptim1Var(),OptimizationOptions(f_abstol = 0.1))
        det_crit = TT(Solvers.x_minimum(res))
        v_crit = TT(Solvers.x_sol(res))
        if det_crit < 0
            detx = det_crit
            vx = v_crit
        else
            return zero(TT)/zero(TT)
        end
    end

    px = is_liquid(phase) ? pmax_sat : pmin_sat

    #the middle det results in an unstable value, that means we can use bisection
    vw = if is_liquid(phase)
        volume(model,px,T,z,phase = :l)::TT
    else
        volume_virial(model,px,T,z)::TT
    end

    xx .= z ./ vw
    detw = det_∂²A∂ϱᵢ²(model,T,xx)
    deta,detb = detw,detx

    if vw < vx
        va,vb = TT(vw),TT(vx)
        deta,detb = TT(detw),TT(detx)
    else
        va,vb = TT(vx),TT(vw)
        deta,detb = TT(detx),TT(detb)
    end

    for i in 1:12
        vc = sqrt(va*vb)
        xx .= z ./ vc
        detc = TT(det_∂²A∂ϱᵢ²(model,T,xx))
        if detc*deta < 0
            va,vb = va,vc
            deta,detb = deta,detc
        elseif detc*detb < 0
            va,vb = vc,vb
            deta,detb = detc,detb
        else
            break
        end
        i >= 7 && detc > 0 && break
    end

    return TT(sqrt(va*vb))
end




"""
    spinodal_temperature(model::EoSModel, p, z; T0, v0, phase)

Calculates the spinodal temperature `T` and volume `V` for a given pressure `p` and composition `z`. Returns a tuple, containing:
- Spinodal temperature `[K]`
- Spinodal volume `[m³]`

Calculates either the liquid or the vapor spinodal point depending on the given starting temperature `T0` and volume `v0` or the `phase`. The keyword `phase` is ignored if `T0` or `v0` is given.
"""
function spinodal_temperature(model::EoSModel,p,z=SA[1.];T0=nothing,v0=nothing,phase=:unknown)
    check_arraysize(model,z)

    is_liquid(phase) || is_vapour(phase) || (!isnothing(v0) && !isnothing(T0)) || error("Either `T0` and `v0` or `phase` have to be specified!")
    model, idx_r = index_reduction(model,z)
    x = z[idx_r]

    # Determine initial guess (if not provided)
    if isnothing(T0) || isnothing(v0)
        v00,T00 = x0_spinodal_temperature(model,p,x,phase)
        _T0 = isnothing(T0) ? T00 : T0*one(T00)
        _v0 = isnothing(v0) ? v00 : v0*one(v00)
    else
        _1 = one(Base.promote_eltype(model,p,z))
        _v0 = _1*v0
        _T0 = _1*T0
        if is_unknown(phase)
            phase = VT_identify_phase(model,_v0,_T0,x)
        end
    end

    # Solve spinodal condition
    function f(X)
        v,T = exp(X[1]),X[2]
        ϱ = x./v
        F1 = det_∂²A∂ϱᵢ²(model, T, ϱ)
        F2 = pressure(model,v,T,x) - p
        return SVector(F1,F2)
    end

    xsol = Solvers.nlsolve2(f,SVector(log(_v0),_T0),Solvers.Newton2Var())
    T_spin = xsol[2]
    v_spin = exp(xsol[1])
    return T_spin,v_spin
end

function x0_spinodal_temperature(model, p, z = SA[1.0], phase = :unknown)

    is_liquid(phase) || is_vapour(phase) || error("`phase` have to be specified!")

    is_vapour(phase) && p < 0 && throw(DomainError(p,"cannot calculate vapour spinodal with negative pressures."))

    TT = Base.promote_eltype(model,p,z)
    nan = zero(TT)/zero(TT)

    if p < 0 && is_liquid(phase)
        Ti = T_scale(model,z)*one(p)
        vi = Ti
        for i in 1:20
            vi = volume(model,zero(p),Ti,z,phase = :l)
            if !isnan(v)
                return vi,Ti
            else
                Ti = 0.9*Ti
            end
        end
        return nan,nan
    end

    v0_edge,dpdT = x0_edge_temperature(model,p,z)
    Tmin_sat,Tmax_sat = extrema(xx -> T_from_dpdT(xx,p),dpdT)
    T0_bubble,T0_dew = v0_edge
    edge,crit,status = _edge_temperature(model,p,z,v0_edge)
    T_edge,v_l,v_v = edge

    TT = Base.promote_eltype(model,p,T_edge,z)
    xx = similar(z,TT)

    status == :failure && (return nan,nan)

    Tx = status == :supercritical ? crit[1] : T_edge

    vx = if status == :supercritical
        volume(model,p,Tx,z,phase = :l, vol0 = crit[3])
    else
        if is_liquid(phase)
            v_l
        else
            v_v
        end
    end

    xx .= z ./ vx

    detx = det_∂²A∂ϱᵢ²(model,Tx,xx)
    status != :supercritical && detx > 0 && (return TT(vx),Tx)

    if status == :supercritical && detx > 0

        v_ref = Ref{TT}(vx)

        f(τ) = let model = model,p = p,z = z,xx  = xx
            tt = 1/τ
            _v = volume(model,p,tt,z,phase = :l,vol0 = v_ref[])
            v_ref[] = _v
            xx .= z ./ _v
            det_∂²A∂ϱᵢ²(model,tt,xx)::TT
        end

        res = Solvers.optimize(f,(1/Tmax_sat,1/Tmin_sat),Solvers.BoundOptim1Var(),OptimizationOptions(f_abstol = 0.1))
        det_crit = TT(Solvers.x_minimum(res))
        τ_crit = TT(Solvers.x_sol(res))
        T_crit = 1/τ_crit
        if det_crit < 0
            detx = det_crit
            Tx = T_crit
            vx = v_ref[]
        else
            return nan,nan
        end
    end

    Tw = is_liquid(phase) ? Tmin_sat : Tmax_sat

    #the middle det results in an unstable value, that means we can use bisection
    vw = if is_liquid(phase)
        volume(model,p,Tw,z,phase = :l)::TT
    else
        volume_virial(model,p,Tw,z)::TT
    end

    xx .= z ./ vw
    detw = det_∂²A∂ϱᵢ²(model,Tw,xx)
    deta,detb = detw,detx

    if vw < vx
        Ta,Tb = TT(Tw),TT(Tx)
        va,vb = TT(vw),TT(vx)
        deta,detb = TT(detw),TT(detx)
    else
        Ta,Tb = TT(Tx),TT(Tw)
        va,vb = TT(vx),TT(vw)
        deta,detb = TT(detx),TT(detb)
    end

    for i in 1:12
        Tc_falsi = Ta - deta*(Tb - Ta)/(detb - deta)
        if Ta <= Tc_falsi <= Tb
            Tc = Tc_falsi
        else
            Tc = 2/(1/Ta + 1/Tb)
        end
        #Tc = 2/(1/Ta + 1/Tb)
        vol0c = is_liquid(phase) ? va : vb
        vc = volume(model,p,Tc,z,phase = phase,vol0 = vol0c)

        xx .= z ./ vc
        detc = TT(det_∂²A∂ϱᵢ²(model,Tc,xx))
        if detc*deta < 0
            Ta,Tb = Ta,Tc
            va,vb = va,vc
            deta,detb = deta,detc
        elseif detc*detb < 0
            Ta,Tb = Tc,Tb
            va,vb = vc,vb
            deta,detb = detc,detb
        else
            break
        end
        abs(deta-detb) < 1e-9 && break
        abs(Ta-Tb) < 1e-9 && break
        i >= 7 && detc > 0 && break
    end

    T = deta < 0 ? Tb : Ta
    v = deta < 0 ? vb : va
    return v,T
end


# Objective function for spinodal calculation -> det(∂²A/∂ϱᵢ) = 0
function det_∂²A∂ϱᵢ²(model,T,ϱ)
    H = Ψ_hessian(model,T,ϱ)
    TT =eltype(H)
    fac = Solvers.unsafe_LU!(H)
    # calculates det(∂²A∂xᵢ² ⋅ ϱ) at V,T constant (see www.doi.org/10.1016/j.fluid.2017.04.009)
    return TT(det(fac))
end

function det_∂²A∂ϱᵢ²(model,V,T,z)
    ϱ = z./V
    return det_∂²A∂ϱᵢ²(model,T,ϱ)
end

liquid_spinodal_zero_limit(model) = liquid_spinodal_zero_limit(model,SA[1.0])
function liquid_spinodal_zero_limit(model::EoSModel,z)
    T = T_scale(model,z)
    p = zero(T)
    Tmax,Tmin = T,T
    nan = zero(T)/zero(T)
    vl = nan
    for i in 1:20
        vl = volume(model,p,T,z,phase =:l)
        if !isnan(vl)
            break
        end
        Tmax = T
        Tmin = 0.9*T
        T = 0.9*T
    end
    if isnan(vl)
        return nan,nan
    end

    for i in 1:30
        T = 0.5*(Tmin + Tmax)
        dT = Tmax - Tmin
        vlx = volume(model,p,T,z,vol0 = vl,phase = :l)
        if isnan(vlx) #T -> Tmax
            Tmax = T
            Tmin = Tmin
            vl = vl
        else #Tmin = T
            Tmin = T
            Tmax = Tmax
            vl = vlx
        end
        dT < 1e-1 && break
    end
    ps = p_scale(model,z)
    lb_v = lb_volume(model,T,z)
    F0(xx) = Obj_crit_pure_sp(xx,model,p,z,ps,lb_v)
    res = Solvers.nlsolve2(F0,SVector(log(vl),T),Solvers.Newton2Var())
    Tw = res[2]
    vw = exp(res[1])
    pw = pressure(model,vw,Tw,z)
    return Tw,vw
end

"""
    spinodal_maximum(model::EoSModel,z;v0=x=x0_crit_mix(model,z))

Returns the maximum temperature/pressure at which there is a mixture spinodal point.

Returns a tuple, containing:
- Spinodal maximum Temperature `[K]`
- Spinodal maximum Pressure `[Pa]`
- Volume at spinodal maximum `[m³]`
"""
function spinodal_maximum(model::EoSModel,z;v0=nothing)
    ∑z = sum(z)
    model_r,idx_r = index_reduction(model,z)

    if length(model_r)==1
        (T_c,p_c,V_c) = crit_pure(model_r)
        return (T_c,p_c,V_c*∑z)
    end

    z_r = z[idx_r]
    z_r ./= ∑z
    if v0 === nothing
        v0 = x0_crit_pure(model_r,z_r)
    end
    Ts = T_scale(model_r,z_r)
    x0 = SVector(v0[1],v0[2]) #could replace for MVector{2}
    f(x) = obj_spinodal_maximum(model_r, exp10(x[2]), Ts*x[1], z_r)
    sol  = Solvers.nlsolve2(f,x0,Solvers.Newton2Var())
    T_c = sol[1]*Ts
    V_c = exp10(sol[2])
    p_c = pressure(model_r, V_c, T_c, z_r)
    return (T_c, p_c, ∑z*V_c)
end

function obj_spinodal_maximum(model,V,T,z)
    f(V) = det_∂²A∂ϱᵢ²(model,V,T,z)
    fv,dfv = Solvers.f∂f(f,V)
    return SVector(fv,dfv)
end

#=
function crit_pure_sp(model)
    T,vl = liquid_spinodal_zero_limit(model)
    vv = zero(vl)
    p = zero(vl)
    ps = p_scale(model,SA[1.0])
    lb_v = lb_volume(model,T,SA[1.0])
    for i in 1:20
        vv = pure_spinodal(model,T,phase = :v)
        p = pressure(model,vv,T)
        F0(w) = Obj_crit_pure_sp(w,model,p,SA[1.0],ps,lb_v)
        res = Solvers.nlsolve2(F0,SVector(log(vl),T),Solvers.Newton2Var())
        T = res[2]
        vl = exp(res[1])
    end
    return T,p,vl,vv
end =#

function Obj_crit_pure_sp(xx,model,p,z,ps,lb_v)
    Vᵢ,Tᵢ = exp(xx[1]),xx[2]
    pᵢ,dpdVᵢ = p∂p∂V(model,Vᵢ,Tᵢ,z)
    return SVector((pᵢ-p)/ps,dpdVᵢ/(ps*lb_v))
end

function eigmin_minimum_pressure(model,T,z,v0hi,v0lo = -second_virial_coefficient(model,T,z))
    f0(v) = diffusive_eigvalue(model,exp(v),T,z)
    ln_vhi = log(v0hi)
    ln_vlo = log(v0lo)
    ln_v = Solvers.optimize(f0,(ln_vhi,ln_vlo),Solvers.BoundOptim1Var())
    v = exp(ln_v)
    return pressure(model,v,T,z),v,f0(ln_v)
end

export spinodal_pressure, spinodal_temperature, spinodal_maximum
