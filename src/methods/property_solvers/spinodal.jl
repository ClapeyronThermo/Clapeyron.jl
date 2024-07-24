
"""
    spinodal_pressure(model::EoSModel, T, x; v0, phase)

Calculates the spinodal pressure and volume for a given temperature and composition. 
Returns a tuple, containing:
- spinodal pressure [`Pa`]
- spinodal volume [`m³`]    
    
Calculates either the liquid or the vapor spinodal point depending on the given starting volume `v0` or the `phase`. The keyword `phase` is ignored if `v0` is given.
"""
function spinodal_pressure(model::EoSModel,T,x=SA[1.];v0=nothing,phase=:unknown,kwargs...)
    x = x/sum(x)
    model, idx_r = index_reduction(model,x)
    x = x[idx_r]
    
    # Determine initial guess (if not provided)
    if isnothing(v0)
        if is_liquid(phase)
            v0 = bubble_pressure(model,T,x;kwargs...)[2]
        elseif is_vapour(phase)
            v0 = dew_pressure(model,T,x;kwargs...)[3]
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

export spinodal_pressure