
"""
    spinodal_pressure(model::EoSModel, T, x; v0, phase)

Calculates the spinodal pressure and volume for a given temperature and composition. Returns a tuple, containing:
- spinodal pressure [`Pa`]
- spinodal volume [`m³`]    
    
Calculates either the liquid or the vapor spinodal point depending on the given starting volume `v0` or the `phase`. The keyword `phase` is ignored if `v0` is given.
"""
function spinodal_pressure(model::EoSModel,T,z=SA[1.];v0=nothing,phase=:unknown)
    TT = Base.promote_eltype(model,T,z)
    V_spin::TT = zero(TT)
    if isnothing(v0) && !is_unknown(phase) && length(model) == 1
        V_spin = pure_spinodal(model,T;phase)
        return pressure(model,V_spin,T),V_spin
    end

    x = z/sum(z)
    model, idx_r = index_reduction(model,x)
    x = x[idx_r]
    
    if isnothing(v0) && !is_unknown(phase) && length(model) == 1
        V_spin = pure_spinodal(model,T;phase)
        return pressure(model,V_spin,T),V_spin
    end
    # Determine initial guess (if not provided)
    if isnothing(v0)
        if is_liquid(phase)
            v0 = bubble_pressure(model,T,x)[2]
        elseif is_vapour(phase)
            v0 = dew_pressure(model,T,x)[3]
        else
            error("Either `v0` or `phase` has to be specified!")
        end
    end

    # Solve spinodal condition
    f!(F,vz) = det_∂²A∂ϱᵢ²(model,F,exp10(vz[1]),T,x)
    r = Solvers.nlsolve(f!,[log10(v0)],LineSearch(Newton()),NEqOptions(), ForwardDiff.Chunk{1}())
    V_spin = exp10(Solvers.x_sol(r)[1])

    if r.info.best_residual[1] < r.options.f_abstol # converged
        return pressure(model,V_spin,T,x), V_spin
    else                                            # not converged
        nan = zero(TT),zero(TT)
        return nan, nan
    end
end

"""
    spinodal_temperature(model::EoSModel, p, x; T0, v0, phase)

Calculates the spinodal pressure and volume for a given pressure and composition. Returns a tuple, containing:
- spinodal temperataure [`K`]
- spinodal volume [`m³`]    

Calculates either the liquid or the vapor spinodal point depending on the given starting temperature `T0` and volume `v0` or the `phase`. The keyword `phase` is ignored if `T0` or `v0` is given.
"""
function spinodal_temperature(model::EoSModel,p,x=SA[1.];T0=nothing,v0=nothing,phase=:unknown)
    x = x/sum(x)
    model, idx_r = index_reduction(model,x)
    x = x[idx_r]
    
    # Determine initial guess (if not provided)
    if isnothing(T0) || isnothing(v0)
        if is_liquid(phase)
            Tv0 = bubble_temperature(model,p,x)[[1,2]]
        elseif is_vapour(phase)
            Tv0 = dew_temperature(model,p,x)[[1,3]]
        else
            error("Either `T0` and `v0` or `phase` have to be specified!")
        end
        T0 = isnothing(T0) ? Tv0[1] : T0
        v0 = isnothing(v0) ? Tv0[2] : v0
    end

    # Solve spinodal condition
    f!(F,Tz) = det_∂²A∂ϱᵢ²(model, F, volume(model, p, Tz[1], x; phase=phase, vol0=v0), Tz[1], x)
    r = Solvers.nlsolve(f!,[T0],LineSearch(Newton()),NEqOptions(;f_abstol=1e-6), ForwardDiff.Chunk{1}())
    T_spin = Solvers.x_sol(r)[1]

    if all(r.info.best_residual .< r.options.f_abstol)  # converged
        return T_spin, volume(model, p, T_spin, x; phase=phase, vol0=v0)
    else                                                # not converged
        return NaN, NaN
    end
end

# Objective function for spinodal calculation -> det(∂²A/∂ϱᵢ) = 0
function det_∂²A∂ϱᵢ²(model,F,v,T,x)
    # calculates det(∂²A∂xᵢ² ⋅ ϱ) at V,T constant (see www.doi.org/10.1016/j.fluid.2017.04.009)
    Av = ϱi -> eos(model, v, T, v.*ϱi)./v
    F[1] = det(ForwardDiff.hessian(Av,x./v))
    return F
end

export spinodal_pressure, spinodal_temperature