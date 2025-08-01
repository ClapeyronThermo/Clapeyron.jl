
"""
    spinodal_pressure(model::EoSModel, T, x; v0, phase)

Calculates the spinodal pressure `p` and `v` volume for a given temperature `T` and composition `z`. Returns a tuple, containing:
- Spinodal pressure `[Pa]`
- Spinodal volume `[m³]`    
    
Calculates either the liquid or the vapor spinodal point depending on the given starting volume `v0` or the `phase`. The keyword `phase` is ignored if `v0` is given.
"""
function spinodal_pressure(model::EoSModel,T,z=SA[1.];v0=nothing,phase=:unknown)
    check_arraysize(model,z)
    ∑z = sum(z)
    if isnothing(v0) && !is_unknown(phase) && length(model) == 1
        V_spin = pure_spinodal(model,T;phase)
        return pressure(model,V_spin,T),V_spin*∑z
    end
    
    xx = z/∑z
    model, idx_r = index_reduction(model,xx)
    x = xx[idx_r]
    
    #reduced model has length = 1
    if isnothing(v0) && !is_unknown(phase) && length(model) == 1
        V_spin = pure_spinodal(model,T;phase)
        return pressure(model,V_spin,T),V_spin*∑z
    end

    # Determine initial guess (if not provided)
    if !isnothing(v0)
        _v0 = one(Base.promote_eltype(model,T,z))*v0
    else
        if is_liquid(phase)
            _v0 = bubble_pressure(model,T,x)[2]
        elseif is_vapour(phase)
            _v0 = dew_pressure(model,T,x)[3]
        else
            error("Either `v0` or `phase` has to be specified!")
        end
    end

    f(vz) = det_∂²A∂ϱᵢ²(model, T, x ./ exp(vz))
    prob = Roots.ZeroProblem(Solvers.to_newton(f),log(_v0))
    log_V_spin = Roots.solve(prob,Roots.Newton())
    V_spin = exp(log_V_spin)
    p_spin = pressure(model,V_spin,T,x)
    return p_spin,V_spin*∑z
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
    ∑z = sum(z)
    xx = z/∑z
    model, idx_r = index_reduction(model,xx)
    x = xx[idx_r]
    
    # Determine initial guess (if not provided)
    if isnothing(T0) || isnothing(v0)
        if is_liquid(phase)
            if p > 0
                T00,v00,_,_ = bubble_temperature(model,p,x)
            else #bubble temperatures are defined for p > 0 
                T00,v0limit = liquid_spinodal_zero_limit(model,x)
                v00  = 0.5*(lb_volume(model,T00,z) + v0limit)
            end
        elseif is_vapour(phase)
            if p <= 0
                throw(DomainError(p,"cannot calculate vapour spinodal with negative pressures."))
            end
            T00,_,v00,_ = dew_temperature(model,p,x)
        else
            error("Either `T0` and `v0` or `phase` have to be specified!")
        end
        _T0 = isnothing(T0) ? T00 : T0*one(T00)
        _v0 = isnothing(v0) ? v00 : v0*one(v00)
    else
        _1 = one(Base.promote_eltype(model,p,z))
        _v0 = _1*v0
        _T0 = _1*T0
    end
    if is_unknown(phase)
        phase = VT_identify_phase(model,_v0,_T0,z)
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

# Objective function for spinodal calculation -> det(∂²A/∂ϱᵢ) = 0
function det_∂²A∂ϱᵢ²(model,T,ϱ)
    # calculates det(∂²A∂xᵢ² ⋅ ϱ) at V,T constant (see www.doi.org/10.1016/j.fluid.2017.04.009)
    Av(ϱi) = Ψ_eos(model, T, ϱi)
    H = ForwardDiff.hessian(Av,ϱ)
    return det(Symmetric(H))
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

export spinodal_pressure, spinodal_temperature
