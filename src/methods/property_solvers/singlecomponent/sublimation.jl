function obj_sublimation_pressure(model::CompositeModel,F,T,vs,vv,p̄,T̄)
    eos_solid(V) = eos(model.solid,V,T,z)
    eos_fluid(V) = eos(model.fluid,V,T,z)
    A_v,Av_v = Solvers.f∂f(eos_fluid,vv)
    A_s,Av_s =Solvers.f∂f(eos_solid,vs)
    μv = muladd(-vv,Av_v,A_v)
    μs = muladd(-vs,Av_s,A_s)
    ps = - Av_s
    pv = - Av_v
    #=
    μs = VT_chemical_potential(model.solid, vs, T)[1]
    μv = VT_chemical_potential(model.fluid, vv, T)[1]
    ps = pressure(model.solid, vs, T)
    pv = pressure(model.fluid, vv, T) =#
    F[1] = (μs - μv)/R̄/T̄
    F[2] = (ps - pv)/p̄
    return F
end

"""
    psub,vs,vv = sublimation_pressure(model::CompositeModel,T;v0=x0_sublimation_pressure(model,T))

Calculates the sublimation pressure of a `CompositeModel` containing a solid and fluid phase EoS, at a specified pressure.
You can pass a tuple of initial values for the volumes `(vs0,vv0)`.

returns:
- Sublimation Pressure [`Pa`]
- Sublimation solid volume at specified temperature [`m³`]
- Sublimation vapour volume at specified temperature [`m³`]
"""

function sublimation_pressure(model::CompositeModel,T;v0=x0_sublimation_pressure(model,T))
    T̄ = T_scale(model.fluid)
    p̄ = p_scale(model.fluid)
    V0 = vec2(log(v0[1]),log(v0[2]),T)
    f!(F,x) = obj_sublimation_pressure(model,F,T,exp10(x[1]),exp10(x[2]),p̄,T̄)
    results = Solvers.nlsolve(f!,V0)
    x = Solvers.x_sol(results)
    vs = exp(x[1])
    vv = exp(x[2])
    return pressure(model.fluid, vv, T), vs, vv
end

function x0_sublimation_pressure(model,T)
    vs0 = lb_volume(model.solid)*6*1.05/π
    vv0 = lb_volume(model.fluid)*6*1.25/π
    return (vs0,vv0)
end