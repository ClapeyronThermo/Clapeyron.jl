function obj_triple_point(model::CompositeModel,F,T,vs,vl,vv,p_scale,T_scale)
    z = SA[1.0]
    fluid = fluid_model(model)
    solid = solid_model(model)
    eos_solid(V) = eos(model.solid,V,T,z)
    eos_fluid(V) = eos(model.fluid,V,T,z)
    A_v,Av_v = Solvers.f∂f(eos_fluid,vv)
    A_l,Av_l = Solvers.f∂f(eos_fluid,vl)
    A_s,Av_s =Solvers.f∂f(eos_solid,vs)
    μv = muladd(-vv,Av_v,A_v)
    μl = muladd(-vl,Av_l,A_l)
    μs = muladd(-vs,Av_s,A_s)
    pv = - Av_v
    pl = - Av_l
    ps = - Av_s
    F[1] = (μs - μv)/R̄/T_scale
    F[2] = (μl - μv)/R̄/T_scale
    F[3] = (ps - pv)/p_scale
    F[4] = (pl - pv)/p_scale
    return F
end

function x0_triple_point(model::CompositeModel,T0 = 0.65*T_scale(fluid_model(model)))
    fluid = fluid_model(model)
    solid = solid_model(model)
    
    #saturation
    p_sat,vl_sat,vv_sat = saturation_pressure(fluid,T0,crit_retry = false)
    H_v = VT_enthalpy(fluid,vv_sat,T0)
    H_l = VT_enthalpy(fluid,vl_sat,T0)
    H_sat = H_v - H_l
    R = Rgas()
    #sublimation:
    vs_sub,vv_sub = x0_sublimation_pressure(model,T0)
    p_sub = pressure(fluid,vv_sub,T0)
    Hs_v = VT_enthalpy(fluid,vv_sub,T0)
    Hs_s = VT_enthalpy(solid,vs_sub,T0)
    H_sub = Hs_v - Hs_s

    #Clausius Clapeyron:
    #=
    log(Ptriple/P0) = - H0/R * (1/Ttriple - 1/T0)
    log(Ptriple) - ln(P0) = -H0/R /Ttriple  + H0/(R*T0)
    log(Ptriple) = H0/(R*T0) - ln(P0) - H0/R /Ttriple 
    ℙ = log(Ptriple), 𝕋 = 1/Ttriple 
    
    ℙ = H0/(R*T0) - ln(P0) - 𝕋*H0/R
    ℙ = A0 - 𝕋*B0

    for sublimation and saturation, we have different definitions of A0,B0. we then solve
    the linear system:
    ℙ = A0_sub - 𝕋*B0_sub
    ℙ = A0_sat - 𝕋*B0_sat
    0 = (A0_sat - A0_sub) - 𝕋*(B0_sat - B0_sub)
    𝕋 = (A0_sat - A0_sub)/(B0_sat - B0_sub)
    =#
    B0_sat = H_sat/R
    B0_sub = H_sub/R
    A0_sat = B0_sat/T0 + log(p_sat)
    A0_sub = B0_sub/T0 + log(p_sub)

    𝕋 = (A0_sat - A0_sub)/(B0_sat - B0_sub)
    ℙ = A0_sat - B0_sat*𝕋
    T = 1/𝕋
    p = exp(ℙ)
    #@show p,T
    vs0 = volume(model,p,T,phase = :s) |> log
    vl0 = volume(model,p,T,phase = :l) |> log
    vv0 = volume(model,p,T,phase = :v) |> log
    return [T/T_scale(fluid_model(model)),vs0,vl0,vv0]
end
"""
    Tt,pt,vs,vl,vv = triple_point(model::CompositeModel;v0 = x0_triple_point(model))

Calculates the triple point of a `CompositeModel` containing solid and fluid phase EoS.

returns:
- Triple point Temperature [`K`]
- Triple point Pressure [`Pa`]
- solid volume at Triple Point [`m³`]
- liquid volume at Triple Point [`m³`]
- vapour volume at Triple Point [`m³`]
"""
function triple_point(model::CompositeModel;v0 = x0_triple_point(model))
    if isnothing(model.solid) || isnothing(model.fluid)
        throw(ArgumentError("triple_point: model must have solid, liquid and gas phases"))
    end

    T̄  = T_scale(model.fluid)
    p̄  = p_scale(model.fluid)

    Ts = T_scale(model.solid)
    f!(F,x) = obj_triple_point(model,F,x[1]*Ts,exp(x[2]),exp(x[3]),exp(x[4]),p̄,T̄)

    results = Solvers.nlsolve(f!,v0,TrustRegion(Newton(),Dogleg()))

    x = Solvers.x_sol(results)

    Tp = x[1]*Ts
    vs = exp(x[2])
    vl = exp(x[3])
    vv = exp(x[4])
    pp = pressure(model.fluid, vv, Tp)
    return (Tp,pp,vs,vl,vv)
end