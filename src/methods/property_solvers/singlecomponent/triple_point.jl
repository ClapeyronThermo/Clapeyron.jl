function obj_triple_point(model::CompositeModel,F,T,vs,vl,vv,p_scale,T_scale)
    μs = VT_chemical_potential(model.solid, vs, T)[1]
    μl = VT_chemical_potential(model.liquid, vl, T)[1]
    μv = VT_chemical_potential(model.gas, vv, T)[1]

    ps = pressure(model.solid, vs, T)
    pl = pressure(model.liquid, vl, T)
    pv = pressure(model.gas, vv, T)

    F[1] = (μs - μl)/R̄/T_scale
    F[2] = (μl - μv)/R̄/T_scale
    F[3] = (ps - pl)/p_scale
    F[4] = (pl - pv)/p_scale
    return F
end

function x0_triple_point(model::CompositeModel)
    T0 = 0.65
    vs0 = log10(lb_volume(model.solid)*1.05*6/π)
    vl0 = log10(lb_volume(model.liquid)*1.25*6/π)
    vv0 = log10(lb_volume(model.gas)*1000)

    return [T0,vs0,vl0,vv0]
end

function triple_point(model::CompositeModel;v0 = x0_triple_point(model))
    if isnothing(model.solid) || isnothing(model.liquid) || isnothing(model.gas)
        throw(ArgumentError("triple_point: model must have solid, liquid and gas phases"))
    end

    T̄  = T_scale(model.liquid)
    p̄  = p_scale(model.liquid)

    f!(F,x) = obj_triple_point(model,F,x[1]*model.solid.params.epsilon[1],exp10(x[2]),exp10(x[3]),exp10(x[4]),p̄,T̄)

    results = Solvers.nlsolve(f!,v0)

    x = Solvers.x_sol(results)

    Tp = x[1]*model.solid.params.epsilon[1]
    vs = exp10(x[2])
    vl = exp10(x[3])
    vv = exp10(x[4])
    pp = pressure(model.gas, vv, Tp)
    return (Tp,pp,vs,vl,vv)
end