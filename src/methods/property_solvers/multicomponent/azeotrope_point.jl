function azeotrope_pressure(model::EoSModel, T; v0 =nothing)
    TYPE = promote_type(eltype(T),eltype([0.5,0.5]))
#     lb_v = lb_volume(model,x)
    ts = T_scales(model)
    pmix = p_scale(model,[0.5,0.5])
    if v0 === nothing
        v0 = x0_bubble_pressure(model,T,[0.5,0.5])
    end
    len = length(v0[1:end-1])
    #xcache = zeros(eltype(x0),len)
    Fcache = zeros(eltype(v0[1:end-1]),len)
    f! = (F,z) -> Obj_az_pressure(model, F, T, exp10(z[1]), exp10(z[2]), z[3:end],z[3:end],ts,pmix)
    r  =Solvers.nlsolve(f!,v0[1:end-1],LineSearch(Newton()))
    sol = Solvers.x_sol(r)
    v_l = exp10(sol[1])
    v_v = exp10(sol[2])
    y = FractionVector(sol[3:end])
    P_sat = pressure(model,v_l,T,y)
    return (P_sat, v_l, v_v, y)
end

function Obj_az_pressure(model::EoSModel, F, T, v_l, v_v, x, y,ts,ps)
    return μp_equality(model::EoSModel, F, T, v_l, v_v, FractionVector(x), FractionVector(y),ts,ps)
end

function azeotrope_temperature(model,p;T0 = nothing)
    f(z) = Obj_azeotrope_temperature(model,z,p)
    if T0 === nothing
        T0 = x0_azeotrope_temperature(model,p)
    end
    fT = Roots.ZeroProblem(f,T0)
    T = Roots.solve(fT)
    p,v_l,v_v,y = azeotrope_pressure(model,T)
    return T,v_l,v_v,y
end

function Obj_azeotrope_temperature(model,T,p)
    p̃,v_l,v_ll,xx = azeotrope_pressure(model,T)
    return p̃-p
end

function x0_azeotrope_temperature(model,p)
    Ti = _sat_Ti(model,p)
    Tmin,Tmax = extrema(Ti)
    return (0.9*Tmin,1.1*Tmax)
end


