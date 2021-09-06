function UCEP_mix(model::EoSModel;v0=nothing)
    if v0 == nothing
        v0 = x0_UCEP_mix(model)
    end  
    ts = T_scales(model,[0.5,0.5])
    pmix = p_scale(model,[0.5,0.5])
    f! = (F,x) -> Obj_UCEP_mix(model, F, x[1], x[2], exp10(x[3]), exp10(x[4]), x[5],ts,pmix)
    r  = Solvers.nlsolve(f!,v0[1:end],LineSearch(Newton()))
    sol = Solvers.x_sol(r)
    x = FractionVector(sol[1])
    y = FractionVector(sol[2])
    V_l = exp10(sol[3])
    V_v = exp10(sol[4])
    T = sol[5]
    p = pressure(model, V_l, T, x)
    return (T, p, V_l, V_v, x, y)
end

function Obj_UCEP_mix(model::EoSModel,F,x,y,V_l,V_v,T,ts,ps)
    y    = FractionVector(y)
    x    = FractionVector(x)
    f(z) = eos(model,V_l,T,z)
    H(z) = ForwardDiff.hessian(f,z)/8.134/T
    L(z) = det(H(z))
    dL(z) = ForwardDiff.gradient(L,z)
    M(z) = [H(z)[1:end-1,:];transpose(dL(z))]
    
    μ_l = VT_chemical_potential(model,V_l,T,x)
    μ_v = VT_chemical_potential(model,V_v,T,y)
    p_l = pressure(model,V_l,T,x)
    p_v = pressure(model,V_v,T,y)
    
    for i in 1:length(x)
        F[i] = (μ_l[i]-μ_v[i])/(R̄*ts[i])
    end
    F[3] = (p_l-p_v)/ps
    
    F[4] = L(x)
    F[5] = det(M(x))
    return F
end

function x0_UCEP_mix(model::EoSModel)
    T0 = T_scale(model,[0.5,0.5])*1.5
    x0 = 0.5
    y0 = 0.75
    v0 = x0_bubble_pressure(model,T0,[x0,1-x0])
    return [x0,y0,v0[1],v0[2],T0]
end