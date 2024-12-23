
"""
    spinodal_pressure(model::EoSModel, T, x; v0, phase)

Calculates the spinodal pressure and volume for a given temperature and composition. Returns a tuple, containing:
- spinodal pressure [`Pa`]
- spinodal volume [`m³`]    
    
Calculates either the liquid or the vapor spinodal point depending on the given starting volume `v0` or the `phase`. The keyword `phase` is ignored if `v0` is given.
"""
function spinodal_pressure(model::EoSModel,T,z=SA[1.];v0=nothing,phase=:unknown)
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

    function fdf(vz)
        fx,dfx = Solvers.f∂f(f,vz)
        return fx,fx/dfx
    end

    prob = Roots.ZeroProblem(fdf,log(_v0))
    log_V_spin = Roots.solve(prob,Roots.Newton(),atol = 1e-6)
    V_spin = exp(log_V_spin)
    p_spin = pressure(model,V_spin,T,x)
    return p_spin,V_spin*∑z
end

"""
    spinodal_temperature(model::EoSModel, p, x; T0, v0, phase)

Calculates the spinodal pressure and volume for a given pressure and composition. Returns a tuple, containing:
- spinodal temperature [`K`]
- spinodal volume [`m³`]    

Calculates either the liquid or the vapor spinodal point depending on the given starting temperature `T0` and volume `v0` or the `phase`. The keyword `phase` is ignored if `T0` or `v0` is given.
"""
function spinodal_temperature(model::EoSModel,p,z=SA[1.];T0=nothing,v0=nothing,phase=:unknown)
    ∑z = sum(z)
    xx = z/∑z
    model, idx_r = index_reduction(model,xx)
    x = xx[idx_r]
    
    # Determine initial guess (if not provided)
    if isnothing(T0) || isnothing(v0)
        if is_liquid(phase)
            T00,v00,_,_ = bubble_temperature(model,p,x)
        elseif is_vapour(phase)
            T00,_,v00,_ = dew_temperature(model,p,x)
        else
            error("Either `T0` and `v0` or `phase` have to be specified!")
        end
        _T0 = isnothing(T0) ? T00 : T0*one(T00)
        _v0 = isnothing(v0) ? v00 : v0*one(v00)
    end

    # Solve spinodal condition
    vcache = [_v0]
    function f(Tz)
        vol0 = vcache[]
        v = volume(model,p,T,x; phase=phase, vol0=vol0)
        ϱ = x./v
        vcache[1] = Solvers.primalval(v)
        det_∂²A∂ϱᵢ²(model, T, ϱ)
    end

    function fdf(Tz)
        fx,dfx = Solvers.f∂f(f,Tz)
        return fx,fx/dfx
    end

    prob = Roots.ZeroProblem(fdf,_T0)
    T_spin = Roots.solve(prob,Roots.Newton(),atol = 1e-6)
    return T_spin, vcache[1]*∑z
end

# Objective function for spinodal calculation -> det(∂²A/∂ϱᵢ) = 0
function det_∂²A∂ϱᵢ²(model,T,ϱ)
    # calculates det(∂²A∂xᵢ² ⋅ ϱ) at V,T constant (see www.doi.org/10.1016/j.fluid.2017.04.009)
    Av(ϱi) = Ψ_eos(model, T, ϱi)
    return det(ForwardDiff.hessian(Av,ϱ))
end

export spinodal_pressure, spinodal_temperature