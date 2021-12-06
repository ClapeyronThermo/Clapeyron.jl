function azeotrope_pressure(model::EoSModel, T; v0 =nothing)
    TYPE = promote_type(eltype(T),eltype([0.5,0.5]))
#     lb_v = lb_volume(model,x)
    ts = T_scales(model,[0.5,0.5])
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
    x   = FractionVector(x)
    y   = FractionVector(y) #julia magic, check misc.jl
    μ_l = VT_chemical_potential(model,v_l,T,x)
    μ_v = VT_chemical_potential(model,v_v,T,y)
    p_l = pressure(model,v_l,T,x)
    p_v = pressure(model,v_v,T,y)
    for i in 1:length(x)
        F[i] = (μ_l[i]-μ_v[i])/(R̄*ts[i])
    end
    F[end] = (p_l-p_v)/ps
    return F
end

function azeotrope_temperature(model,p)
    f(z) = Obj_azeotrope_temperature(model,z,p)
    pure = split_model(model)
    sat = saturation_temperature.(pure,p)
    Ti   = zeros(length(pure),1)
    for i ∈ 1:length(pure)
        if isnan(sat[i][1])
            Tc,pc,vc = crit_pure(pure[i])
            g(x) = p-pressure(pure[i],vc,x,[1.])
            Ti[i] = Roots.find_zero(g,(Tc))
        else
            Ti[i] = sat[i][1]
        end
    end
    T = Roots.find_zero(f,(minimum(Ti)*0.9,maximum(Ti)*1.1))
    p,v_l,v_v,y = azeotrope_pressure(model,T)
    return T,v_l,v_v,y
end

function Obj_azeotrope_temperature(model,T,p)
    p̃,v_l,v_ll,xx = azeotrope_pressure(model,T)
    return p̃-p
end

