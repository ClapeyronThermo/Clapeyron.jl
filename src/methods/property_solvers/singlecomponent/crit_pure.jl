"""
    crit_pure(model::EoSModel,x0=nothing)
Calculates the critical point of a single component modelled by `model`. 
Returns `(Tc, pc, Vc)` where `Tc` is the critical temperature (in K), `pc` is the critical pressure (in Pa) and `Vc` is the critical volume (in  m³)
"""
function crit_pure(model::EoSModel,x0=nothing;options = NEqOptions())
    single_component_check(crit_pure,model)

    #f! = (F,x) -> obj_crit(model, F, x[1]*T̄, exp10(x[2]))
    if x0 === nothing
        x0 = x0_crit_pure(model)
    end
    x01,x02 = x0
    T̄  = T_scale(model)*one(x01*one(x02))
    x0 = vec2(x01,x02,T̄)
    f! = ObjCritPure(model,T̄,x0)
    solver_res = Solvers.nlsolve(f!, x0, TrustRegion(Newton(approx = NLSolvers.Inverse()), NWI()), options)
    r  = Solvers.x_sol(solver_res)
    T_c = r[1]*T̄
    V_c = exp10(r[2])
    p_c = pressure(model, V_c, T_c)
    return (T_c, p_c, V_c)
end

function ObjCritPure(model,T̄,x0)
    xcache = copy(x0)

    function f!(F,x)
        T_c = x[1]*T̄
        V_c = exp10(x[2])
        ∂²A∂V², ∂³A∂V³ = ∂²³f(model, V_c, T_c, SA[1.0])
        F[1] = -∂²A∂V²
        F[2] = -∂³A∂V³
        return F
    end

    function fj!(F,J,x)
        ForwardDiff.jacobian!(J,f!,F,x)
        return F,J
    end

    function j!(J,x)
        ForwardDiff.jacobian!(J,f!,xcache,x)
        return J
    end

    function jv!(inc)
        return nothing
    end

    return Solvers.NLSolvers.VectorObjective(f!,j!,fj!,jv!) |> Solvers.NLSolvers.NEqProblem
end