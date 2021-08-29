## VLLE Solver
function Obj_three_phase(model::EoSModel, F, T, v_l, v_ll, v_v, x, xx, y,ts,ps)
    x   = FractionVector(x)
    y   = FractionVector(y)
    xx  = FractionVector(xx)#julia magic, check misc.jl
    μ_l   = VT_chemical_potential(model,v_l,T,x)
    μ_ll  = VT_chemical_potential(model,v_ll,T,xx)
    μ_v   = VT_chemical_potential(model,v_v,T,y)
    p_l   = pressure(model,v_l,T,x)
    p_ll  = pressure(model,v_ll,T,xx)
    p_v   = pressure(model,v_v,T,y)
    n_c   = length(model.components)
    for i in 1:n_c
        F[i] = (μ_l[i]-μ_v[i])/(R̄*ts[i])
        F[i+n_c] = (μ_ll[i]-μ_v[i])/(R̄*ts[i])
    end
    F[end-1] = (p_l-p_v)/ps
    F[end]   = (p_ll-p_v)/ps
    return F
end

function three_phase(model::EoSModel, T; v0 =nothing)
#     TYPE = promote_type(eltype(T),eltype(x))
#     lb_v = lb_volume(model,x)
    if v0 === nothing
        v0 = x0_three_phase(model,T)
    end
    ts = T_scales(model,[v0[4],1-v0[4]])
    pmix = p_scale(model,[v0[4],1-v0[4]])
    len = length(v0[1:end])
    #xcache = zeros(eltype(x0),len)
    Fcache = zeros(eltype(v0[1:end]),len)
    f! = (F,z) -> Obj_three_phase(model, F, T, exp10(z[1]), exp10(z[2]), exp10(z[3]), z[4], z[5], z[6],ts,pmix)
    r  = Solvers.nlsolve(f!,v0[1:end],LineSearch(Newton()))
    sol = Solvers.x_sol(r)
    v_l = exp10(sol[1])
    v_ll = exp10(sol[2])
    v_v = exp10(sol[3])
    x   = sol[4]
    xx = sol[5]
    y  = sol[6]
    x = FractionVector(x)
    xx = FractionVector(xx)
    y = FractionVector(y)
    P_sat = pressure(model,v_l,T,x)
    return (P_sat, v_l, v_ll, v_v, x, xx, y)
end

function x0_three_phase(model::EoSModel, T)
    pure = split_model(model)
    sat  = sat_pure.(pure,T)
    x0    = [0.75,0.25]
    xx0   = [0.25,0.75]
    y0    = [0.5,0.5]
    v_l0  = sat[1][2]*x0[1]+sat[2][2]*x0[2]
    v_ll0 = sat[1][2]*xx0[1]+sat[2][2]*xx0[2]
    v_v0  = sat[1][3]*y0[1]+sat[2][3]*y0[2]
    return [log10(v_l0),log10(v_ll0),log10(v_v0),x0[1],xx0[1],y0[1]]
end
