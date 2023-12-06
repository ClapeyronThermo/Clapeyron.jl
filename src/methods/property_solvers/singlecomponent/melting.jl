function obj_melting_pressure(model::CompositeModel,F,T,vs,vl,p̄,T̄)
    μs = VT_chemical_potential(model.solid, vs, T)[1]
    μl = VT_chemical_potential(model.liquid, vl, T)[1]

    ps = pressure(model.solid, vs, T)
    pl = pressure(model.liquid, vl, T)

    F[1] = (μs - μl)/R̄/T̄
    F[2] = (ps - pl)/p̄
    return F
end

function melting_pressure(model::CompositeModel,T;v0=x0_melting_pressure(model,T))
    T̄ = T_scale(model.liquid)
    p̄ = p_scale(model.liquid)

    f!(F,x) = obj_melting_pressure(model,F,T,exp10(x[1]),exp10(x[2]),p̄,T̄)

    results = Solvers.nlsolve(f!,v0)

    x = Solvers.x_sol(results)

    vs = exp10(x[1])
    vl = exp10(x[2])
    return pressure(model.liquid, vl, T), vs, vl
end

function x0_melting_pressure(model,T)
    vs0 = log10(lb_volume(model.solid)*6*1.05/π)
    vl0 = log10(lb_volume(model.liquid)*6*1.25/π)
    return [vs0,vl0]
end