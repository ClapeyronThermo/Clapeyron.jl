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
    #if type !== nothing
    #    T̄ = T̄*oneunit(type)
    #end
    
    _1 = oneunit(T̄)
    #x0 = SVector(_1*x01,_1*x02)
    x0 = vec2(x01,x02,_1)
    f!(F,x) = ObjCritPure(model,F,T̄,x)
    solver_res = Solvers.nlsolve(f!, x0, TrustRegion(Newton(), NLSolvers.NWI()), options)
    r  = Solvers.x_sol(solver_res)
    T_c = r[1]*T̄
    V_c = exp10(r[2])
    p_c = pressure(model, V_c, T_c)
    return (T_c, p_c, V_c)
end

function ObjCritPure(model::T,F,T̄,x) where T
    T_c = x[1]*T̄
    V_c = exp10(x[2])
    ∂²A∂V², ∂³A∂V³ = ∂²³f(model, V_c, T_c, SA[1.0])
    F1 = -∂²A∂V²
    F2 = -∂³A∂V³
    #return SVector(F1,F2)
    F[1] = F1;F[2] = F2; return F
end