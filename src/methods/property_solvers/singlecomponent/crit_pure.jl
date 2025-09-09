"""
    crit_pure(model::EoSModel,x0=nothing)

Calculates the critical point of a single component modelled by `model`. 
Returns a tuple, containing:

- Critical temperature `[K]`
- Critical pressure `[Pa]`
- Critical molar volume `[m³·mol⁻¹]`
"""

function crit_pure(model::EoSModel)
    satmodel = saturation_model(model)
    if satmodel !== model
        return crit_pure(satmodel)
    else
        return crit_pure(model,nothing)
    end
end

function crit_pure(model::EoSModel,x0,z = SA[1.0];options = NEqOptions())
    check_arraysize(model,z)
    #f! = (F,x) -> obj_crit(model, F, x[1]*T̄, exp10(x[2]))
    if x0 === nothing
        x0 = x0_crit_pure(model,z)
    end
    x01,x02 = x0
    T̄  = T_scale(model,z)*one(x01*one(x02))
    #if type !== nothing
    #    T̄ = T̄*oneunit(type)
    #end
    _1 = oneunit(T̄)
    #x0 = SVector(_1*x01,_1*x02)
    x0 = vec2(primalval(x01),primalval(x02*log(10)),primalval(_1))
    primalmodel = primalval(model)
    zz = z/sum(z)
    f!(F,x) = ObjCritPure(primalmodel,F,primalval(T̄),x,zz)
    solver_res = Solvers.nlsolve(f!, x0, TrustRegion(Newton(), NLSolvers.NWI()), options)
    #display(solver_res)
    r  = Solvers.x_sol(solver_res)
    T_c = r[1]*T̄
    V_c = exp(r[2])
    p_c = pressure(model, V_c, T_c, zz)
    crit = (T_c, p_c, V_c)
    return crit_pure_ad(model,crit,z)
end

function crit_pure_ad(model,crit,z)
    if has_dual(model)
        T_c_primal, p_c_primal, V_c_primal = crit
        T̄  = T_scale(model)
        x = SVector(T_c_primal/T̄,log(V_c_primal))
        f(zz) = __ObjCritPure(model,T̄,z,zz)
        F,J = Solvers.J2(f,x)
        ∂x = J\F
        r = x .- ∂x
        T_c = r[1]*T̄
        V_c = exp(r[2])
        P_c = pressure(model,V_c,T_c,z)
        return (T_c,P_c,V_c)
    else
        return crit
    end
end

function ObjCritPure(model::T,F,T̄,x,z) where T
    sv = __ObjCritPure(model,T̄,x,z)
    F[1] = sv[1]
    F[2] = sv[2]
    return F
end

function __ObjCritPure(model::T,T̄,x,z) where T
    T_c = x[1]*T̄
    V_c = exp(x[2])
    RT = Rgas(model)*T_c
    ∂²A∂V²_scale = V_c*V_c/RT
    ∂³A∂V³_scale = ∂²A∂V²_scale*V_c
    ∂²A∂V², ∂³A∂V³ = ∂²³f(model, V_c, T_c, z)
    F1 = -∂²A∂V²*∂²A∂V²_scale
    F2 = -∂³A∂V³*∂³A∂V³_scale
    return SVector((F1,F2))
end

"""
    mechanical_critical_point(model,z = SA[1.0],x0 = x0_crit_pure(model,z))

Returns the mechanical critical point of an `EoSModel`. The mechanical critical point is the point where the first and second derivatives of the pressure function vanish.
For single component models, the function is equivalent to [`crit_pure`](@ref).
Returns a tuple, containing:

- Mechanical critical temperature `[K]`
- Mechanical critical pressure `[Pa]`
- Mechanical critical molar volume `[m³·mol⁻¹]`
"""
function mechanical_critical_point end

mechanical_critical_point(model,z,x0) = crit_pure(model,x0,z)
mechanical_critical_point(model,z) = mechanical_critical_point(model,z,x0_crit_pure(model,z))
mechanical_critical_point(model) = crit_pure(model)