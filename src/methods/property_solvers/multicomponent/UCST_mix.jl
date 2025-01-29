"""
    UCST_mix(model::EoSModel,T;v0=x0_UCST_mix(model,T))

Calculates the Upper critical solution point of a mixture at a given Temperature.

returns:
- UCST Pressure [`Pa`]
- volume at UCST Point [`m³`]
- molar composition at UCST Point
"""
function UCST_mix(model::EoSModel,T;v0=nothing)
    if v0 === nothing
        v0 = x0_UCST_mix(model,T)
    end
    x0 = vcat(v0[1],v0[2][1:end-1])
    f! = (F,x) -> Obj_UCST_mix(model, F, x[2:end], exp10(x[1]), T)
    r  = Solvers.nlsolve(f!,x0,LineSearch(Newton2(x0)))
    sol = Solvers.x_sol(r)
    !all(<(r.options.f_abstol),r.info.best_residual) && (sol .= NaN)
    z_c = FractionVector(sol[2:end])
    V_c = exp10(sol[1])
    p_c = pressure(model, V_c, T, z_c)
    return (p_c, V_c, z_c)
end

"""
    x0_UCST_mix(model::EoSModel,T)

Initial point for `UCST_mix(model,T)`.

Returns a tuple, containing:
- Base 10 logarithm initial guess for liquid composition `[m³]`
- Initial guess for molar fractions at UCST Point (default: equimolar)
"""
function x0_UCST_mix(model::EoSModel,T)
    x0 = Fractions.zeros(length(model))
    p = p_scale(model,x0)
    V  = x0_volume(model,p,T,x0,phase = :l)
    return (log10(V),x0)
end

function Obj_UCST_mix(model::EoSModel,F,z,V,T)
    z    = FractionVector(z)
    L,detM = mixture_critical_constraint(model,V,T,z)
    F[1] = L
    F[2] = detM
    return F
end