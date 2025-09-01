function obj_triple_point(model::CompositeModel,F,T,vs,vl,vv,p_scale,T_scale)
    fluid = fluid_model(model)
    solid = solid_model(model)
    z = SA[1.0]
    RT = Rgas(fluid)*T
    f_solid(V) = a_res(solid,V,T,z)
    f_fluid(V) = a_res(fluid,V,T,z)
    As,Avs = Solvers.f∂f(f_solid,vs)
    Al,Avl = Solvers.f∂f(f_fluid,vl)
    Av,Avv = Solvers.f∂f(f_fluid,vv)
    ps,pl,pv = RT*(-Avs + 1/vs),RT*(-Avl + 1/vl),RT*(-Avv + 1/vv)
    F[1] = As - vs*Avs - Av + vv*Avv + log(vv/vs)
    F[2] = Al - vl*Avl - Av + vv*Avv + log(vv/vl)
    F[3] = (ps - pv)/p_scale
    F[4] = (pl - pv)/p_scale
    return F
end

function obj_triple_point_gibbs(model::CompositeModel,F,T,p,T_scale,k1,k2)
    fluid = fluid_model(model)
    solid = solid_model(model)
    gs = gibbs_energy(solid,p,T,SA[1.0],phase = :solid) + k1 + k2*T
    gl = gibbs_energy(fluid,p,T,SA[1.0],phase = :liquid)
    gv = gibbs_energy(fluid,p,T,SA[1.0],phase = :gas)
    RT = Rgas()*T_scale
    F[1] = (gs - gv)/RT
    F[2] = (gs - gl)/RT
    return F
end

function x0_triple_point(model::CompositeModel,T0 = 0.65*T_scale(fluid_model(model)))
    fluid = fluid_model(model)
    solid = solid_model(model)

    #saturation
    p_sat,vl_sat,vv_sat = saturation_pressure(fluid,T0,crit_retry = false)
    Ksat = -dpdT_saturation(fluid,vl_sat,vv_sat,T0)*T0*T0/p_sat
    vs_sub,vv_sub,p_sub = x0_sublimation_pressure(model,T0)
    if fluid isa GibbsBasedModel || solid isa GibbsBasedModel
        Ksub = -dpdT_saturation_gibbs(solid,fluid,p_sub,T0,phase1 = :solid,phase2 = :liquid)*T0*T0/p_sub
    else
        Ksub = -dpdT_saturation(solid,fluid,vs_sub,vv_sub,T0)*T0*T0/p_sub
    end
    #=

    Clausius Clapeyron
    log(Ptriple/P_sat) = K_sat * (1/Ttriple - 1/Tsat)
    log(Ptriple/P_sub) = K_sub * (1/Ttriple - 1/Tsub)
    Ksat = [dpdT*T*T/p](p = p0,T = T0)
    Tsat = Tsub = T0
    log(Ptriple/P_sat) = K_sat * (1/Ttriple - 1/T0)
    log(Ptriple/P_sub) = K_sub * (1/Ttriple - 1/T0)

    p̃ = log(P), T̃ = 1/T
    we solve the following system of equations:
    p̃(triple) - K_sat*T̃(triple) = -K_sat*T̃0 + p̃(sat)
    p̃(triple) - K_sub*T̃(triple) = -K_sub*T̃0 + p̃(sub)
    =#
    _1 = one(Ksub)
    T̃0 = 1/T0
    p̃sat,p̃sub = log(p_sat),log(p_sub)
    M = @SMatrix [_1 -Ksat
                  _1 -Ksub]
    b = SVector(-Ksat*T̃0 + p̃sat,-Ksub*T̃0 + p̃sub)
    p̃,T̃ = M\b
    p,T = exp(p̃),1/T̃
    vs0 = volume(model,p,T,phase = :s)
    vl0 = volume(model,p,T,phase = :l)
    vv0 = volume(model,p,T,phase = :v)
    return [p,T,vs0,vl0,vv0]
end
"""
    Tt,pt,vs,vl,vv = triple_point(model::CompositeModel;v0 = x0_triple_point(model))

Calculates the triple point of a `CompositeModel` containing solid and fluid phase EoS.

returns:
- Triple point Temperature `[K]`
- Triple point Pressure `[Pa]`
- Solid volume at Triple Point `[m³]`
- Liquid volume at Triple Point `[m³]`
- Vapour volume at Triple Point `[m³]`
"""
function triple_point(model::CompositeModel;v0 = x0_triple_point(model))
    if isnothing(model.solid) || isnothing(model.fluid)
        throw(ArgumentError("triple_point: model must have solid, liquid and gas phases"))
    end
    if model.solid isa GibbsBasedModel || model.fluid isa GibbsBasedModel #the second condition is unlikely.
        return triple_point_gibbs(model,v0)
    else
        return triple_point_helmholtz(model,v0)
    end
    return (Tp,pp,vs,vl,vv)
end

function triple_point_helmholtz(model,v0)
    T̄  = T_scale(model.fluid)
    p̄  = p_scale(model.fluid)
    Ts = T_scale(model.solid)
    x0 = [v0[2]/T̄,log(v0[3]),log(v0[4]),log(v0[5])]
    f!(F,x) = obj_triple_point(model,F,x[1]*Ts,exp(x[2]),exp(x[3]),exp(x[4]),p̄,T̄)
    results = Solvers.nlsolve(f!,x0,TrustRegion(Newton(),Dogleg()))
    x = Solvers.x_sol(results)
    Tp = x[1]*Ts
    vs = exp(x[2])
    vl = exp(x[3])
    vv = exp(x[4])
    pp = pressure(model.fluid, vv, Tp)
    return (Tp,pp,vs,vl,vv)
end

function triple_point_gibbs(model,v0)
    T̄  = T_scale(model.solid)
    x0 = [log(v0[1]),v0[2]/T̄]
    k1,k2 = calculate_gibbs_reference_state(model)
    f!(F,x) = obj_triple_point_gibbs(model::CompositeModel,F,x[2]*T̄,exp(x[1]),T̄,k1,k2)
    results = Solvers.nlsolve(f!,x0,TrustRegion(Newton(),Dogleg()))
    logp,Tr = Solvers.x_sol(results)
    p,T = exp(logp),Tr*T̄
    vs = volume(model.solid,p,T,phase = :solid)
    vl = volume(model.fluid,p,T,phase = :liquid)
    vv = volume(model.fluid,p,T,phase = :vapour)
    return (T,p,vs,vl,vv)
end