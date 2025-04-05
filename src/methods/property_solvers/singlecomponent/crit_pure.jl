"""
    crit_pure(model::EoSModel,x0=nothing)
Calculates the critical point of a single component modelled by `model`. 
Returns `(Tc, pc, Vc)` where `Tc` is the critical temperature (in K), `pc` is the critical pressure (in Pa) and `Vc` is the critical volume (in  m³)
"""

function crit_pure(model::EoSModel)
    satmodel = saturation_model(model)
    if satmodel !== model
        return crit_pure(satmodel)
    else
        return crit_pure(model,nothing)
    end
end

function crit_pure(model::EoSModel,x0;options = NEqOptions())
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
    x0 = vec2(primalval(x01),primalval(x02*log(10)),primalval(_1))
    primalmodel = primalval(model)
    f!(F,x) = ObjCritPure(primalmodel,F,primalval(T̄),x)
    solver_res = Solvers.nlsolve(f!, x0, TrustRegion(Newton(), NLSolvers.NWI()), options)
    #display(solver_res)
    r  = Solvers.x_sol(solver_res)
    T_c = r[1]*T̄
    V_c = exp(r[2])
    p_c = pressure(model, V_c, T_c)

    crit = (T_c, p_c, V_c)
    return crit_pure_ad(model,crit)
end

function crit_pure_ad(model,crit)
    if has_dual(model)
        T_c_primal, p_c_primal, V_c_primal = crit
        T̄  = T_scale(model)
        x = SVector(T_c_primal/T̄,log(V_c_primal))
        f(z) = __ObjCritPure(model,T̄,z)
        F,J = Solvers.J2(f,x)
        ∂x = J\F
        r = x .- ∂x
        T_c = r[1]*T̄
        V_c = exp(r[2])
        P_c = pressure(model,V_c,T_c)
        return (T_c,P_c,V_c)
    else
        return crit
    end
end

function ObjCritPure(model::T,F,T̄,x) where T
    sv = __ObjCritPure(model,T̄,x)
    F[1] = sv[1]
    F[2] = sv[2]
    return F
end

function __ObjCritPure(model::T,T̄,x) where T
    T_c = x[1]*T̄
    V_c = exp(x[2])
    RT = Rgas(model)*T_c
    ∂²A∂V²_scale = V_c*V_c/RT
    ∂³A∂V³_scale = ∂²A∂V²_scale*V_c
    ∂²A∂V², ∂³A∂V³ = ∂²³f(model, V_c, T_c, SA[1.0])
    F1 = -∂²A∂V²*∂²A∂V²_scale
    F2 = -∂³A∂V³*∂³A∂V³_scale
    return SVector((F1,F2))
end