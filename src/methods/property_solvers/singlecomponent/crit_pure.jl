"""
    crit_pure(model::EoSModel,x0=nothing)

Calculates the critical point of a single component modelled by `model`. 

Returns `(Tc, pc, Vc)` where `Tc` is the critical temperature (in K), `pc` is the critical pressure (in Pa) and `Vc` is the critical volume (in  m³)
"""
function crit_pure(model::EoSModel,x0=nothing)
    !isone(length(model)) && throw(error("$model have more than one component."))
    T̄  = T_scale(model)
    f! = (F,x) -> obj_crit(model, F, x[1]*T̄, exp10(x[2]))
    if x0 === nothing
        x0 = x0_crit_pure(model)
    end
    x0 = MVector((x0[1],x0[2]))
    #x0 = [x0[1],x0[2]]
    solver_res = Solvers.nlsolve(f!, x0)
    #print(solver_res)
    r  = Solvers.x_sol(solver_res)
    T_c = r[1]*T̄
    V_c = exp10(r[2])
    p_c = pressure(model, V_c, T_c)
    return (T_c, p_c, V_c)
end

function obj_crit(model::EoSModel, F, T_c, V_c)
    _1 = one(T_c+V_c)
    ∂²A∂V², ∂³A∂V³ = ∂²³f(model, V_c, T_c, SA[_1])
    F[1] = -∂²A∂V²
    F[2] = -∂³A∂V³
    return F
end
