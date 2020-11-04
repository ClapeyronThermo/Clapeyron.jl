## Standard pressure solver
function get_volume(model::EoS, p, T, z=[1.]; phase = "unknown")
    z = create_z(model, z)
    components = model.components

    N = length(p)

    ub = [Inf]
    lb = lb_volume(model,z; phase = phase)

    x0 = x0_volume(model,z; phase = phase)

    Vol = []
    if phase == "unknown"
        for i in 1:N
            f(v) = eos(model, z, 10^v[1], T[i]) + 10^v[1]*p[i]
            (f_best,v_best) = Solvers.tunneling(f,lb,ub,x0)
            append!(Vol,10^v_best[1])
        end
    else
        opt_min = NLopt.Opt(:LD_MMA, length(ub))
        opt_min.lower_bounds = lb
        opt_min.upper_bounds = ub
        opt_min.xtol_rel     = 1e-8
        for i in 1:N
            f(v)   = eos(model, z, 10^v[1], T[i]) + 10^v[1]*p[i]
            obj_f0 = x -> f(x)
            obj_f  = (x,g) -> Solvers.NLopt_obj(obj_f0,x,g)
            opt_min.min_objective =  obj_f
            (f_min,v_min) = NLopt.optimize(opt_min, x0)
            append!(Vol, 10^v_min[1])
        end
    end
    return Vol
end


## Pure saturation conditions solver
function get_sat_pure(model::EoS, T)
    components = model.components
    v0    = x0_sat_pure(model)
    v_l   = []
    v_v   = []
    P_sat = []
    for i in 1:length(T)
        f! = (F,x) -> Obj_Sat(model, F, T[i], 10^x[1], 10^x[2])
        j! = (J,x) -> Jac_Sat(model, J, T[i], 10^x[1], 10^x[2])
        r  =nlsolve(f!,j!,v0)
        append!(v_l,10^r.zero[1])
        append!(v_v,10^r.zero[2])
        append!(P_sat,get_pressure(model,v_v[i],T[i]))
        v0 = r.zero
    end
    return (P_sat, v_l, v_v)
end

function Obj_Sat(model::EoS, F, T, v_l, v_v)
    components = model.components
    fun(x) = eos(model, create_z(model, [x[1]]), x[2], T)
    df(x)  = ForwardDiff.gradient(fun,x)
    df_l = df([1,v_l[1]])
    df_v = df([1,v_v[1]])
    (p_scale,μ_scale) = scale_sat_pure(model)
    F[1] = (df_l[2]-df_v[2])*p_scale
    F[2] = (df_l[1]-df_v[1])*μ_scale
end

function Jac_Sat(model::EoS, J, T, v_l, v_v)
    components = model.components
    fun(x) = eos(model, create_z(model, [x[1]]), x[2], T)
    d2f(x) = ForwardDiff.hessian(fun,x)
    d2f_l = d2f([1,v_l[1]])
    d2f_v = d2f([1,v_v[1]])
    (p_scale,μ_scale) = scale_sat_pure(model)
    J[1,1] =  v_l[1]*d2f_l[2,2]*log(10)*p_scale
    J[1,2] = -v_v[1]*d2f_v[2,2]*log(10)*p_scale
    J[2,1] =  v_l[1]*d2f_l[1,2]*log(10)*μ_scale
    J[2,2] = -v_v[1]*d2f_v[1,2]*log(10)*μ_scale
end

function get_enthalpy_vap(model::EoS, T)
    (P_sat,v_l,v_v) = get_sat_pure(model,T)
    fun(x) = eos(model, create_z(model,[1.0]), x[1], x[2])
    df(x)  = ForwardDiff.gradient(fun,x)
    H_vap = []
    for i in 1:length(T)
        H_l = fun([v_l[i],T[i]])-df([v_l[i],T[i]])[2]*T[i]-df([v_l[i],T[i]])[1]*v_l[i]
        H_v = fun([v_v[i],T[i]])-df([v_v[i],T[i]])[2]*T[i]-df([v_v[i],T[i]])[1]*v_v[i]
        append!(H_vap,H_v-H_l)
    end
    return H_vap
end
## Pure critical point solver
function get_crit_pure(model::EoS; units = false, output=[u"K", u"Pa", u"m^3"])
    components = model.components
    T̄  = T_crit_pure(model)
    f! = (F,x) -> Obj_Crit(model, F, x[1]*T̄, 10^x[2])
    # j! = (J,x) -> Jac_Crit(J,eos,model,x[1]*model.params.epsilon[(1, 1)],10^x[2])
    x0 = x0_crit_pure(model)
    r  = nlsolve(f!,x0)
    T_c = r.zero[1]*T̄
    v_c = 10^r.zero[2]
    p_c = get_pressure(model, v_c, T_c)
    if units
        return (uconvert(output[1], T_c*u"K"), uconvert(output[2], p_c*u"Pa"), uconvert(output[2], v_c*u"m^3"))
    else
        return (T_c, p_c, v_c)
    end
end

function Obj_Crit(model::EoS, F, T_c, v_c)
    fun(x)  = eos(model, create_z(model, [1]), x[1], T_c)
    df(x)   = ForwardDiff.derivative(fun,x)
    d2f(x)  = ForwardDiff.derivative(df,x)
    d3f(x)  = ForwardDiff.derivative(d2f,x)
    F[1] = d2f(v_c)
    F[2] = d3f(v_c)
end

## Mixture saturation solver
function get_sat_mix_Tx(model, T, x)
    components = model.components
    y0    = 10 .*x[1,:]./(1 .+x[1,:].*(10-1))
    y0    = y0 ./sum(y0[i] for i in 1:length(components))
    X     = create_z(model,x[1,:])
    Y0    = create_z(model,y0)
    v0    = [log10(π/6*N_A*sum(X[i]*model.params.segment[i]*model.params.sigma[i]^3 for i in components)/0.45),
             log10(π/6*N_A*sum(Y0[i]*model.params.segment[i]*model.params.sigma[i]^3 for i in components)/1e-4)]
    append!(v0,y0[1:end-1])
    v_l   = []
    v_v   = []
    y     = deepcopy(x)
    P_sat = []
    for i in 1:size(x)[1]
        f! = (F,z) -> Obj_Sat_mix_Tx(model, F, T, 10^z[1], 10^z[2], x[i,:], z[3:end])
        j! = (J,z) -> Jac_Sat_mix_Tx(model, J, T, 10^z[1], 10^z[2], x[i,:], z[3:end])
        r  =nlsolve(f!,j!,v0)
        append!(v_l,10^r.zero[1])
        append!(v_v,10^r.zero[2])
        append!(P_sat,get_pressure(model,v_l[i],T,x[i,:]))
        y[i,1:end-1] = r.zero[3:end]
        y[i,end] = 1-sum(r.zero[3:end])
        v0 = r.zero
    end
    return (P_sat, v_l, v_v, y)
end

function Obj_Sat_mix_Tx(model, F, T, v_l, v_v, x, y)
    components = model.components
    append!(y,1-sum(y[i] for i in 1:(length(components)-1)))

    fun(z) = eos(model, create_z(model, z[1:end-1]), z[end], T)
    df(z)  = ForwardDiff.gradient(fun,z)
    X = deepcopy(x)
    Y = deepcopy(y)
    df_l = df(append!(X,v_l))
    df_v = df(append!(Y,v_v))
    for i in 1:length(components)
        F[i] = (df_l[i]-df_v[i])/R̄/model.params.epsilon[components[i]]
    end
    F[end] = (df_l[end]-df_v[end])*model.params.sigma[components[1]]^3*N_A/R̄/model.params.epsilon[components[1]]
end

function Jac_Sat_mix_Tx(model, J, T, v_l, v_v, x, y)
    components = model.components
    append!(y,1-sum(y[i] for i in 1:(length(components)-1)))

    fun(z) = eos(model, create_z(model, z[1:end-1]), z[end], T)
    df(z)  = ForwardDiff.gradient(fun,z)
    d2f(z) = ForwardDiff.hessian(fun,z)
    X = deepcopy(x)
    Y = deepcopy(y)
    d2f_l = d2f(append!(X,v_l))
    d2f_v = d2f(append!(Y,v_v))
    for i in 1:(length(components))
        J[i,1] =  log(10)*v_l*d2f_l[i,end]/R̄/model.params.epsilon[components[i]]
        J[i,2] = -log(10)*v_v*d2f_v[i,end]/R̄/model.params.epsilon[components[i]]
    end
    J[end,1] = log(10)*v_l*d2f_l[end,end]*model.params.sigma[components[1]]^3*N_A/R̄/model.params.epsilon[components[1]]
    J[end,2] =-log(10)*v_v*d2f_v[end,end]*model.params.sigma[components[1]]^3*N_A/R̄/model.params.epsilon[components[1]]

    for j in 1:(length(components)-1)
        J[end,j+2] = (d2f_v[end,end-1]-d2f_v[end,j])*model.params.sigma[components[1]]^3*N_A/R̄/model.params.epsilon[components[1]]
    end

    for i in 1:(length(components))
        for j in 1:(length(components)-1)
                J[i,j+2]= -(d2f_v[i,j]-d2f_v[i,end-1])/R̄/model.params.epsilon[components[i]]
        end
    end
end


## Mixture critical point solver
# function get_crit_mix(model::SAFT,x_c)
#     components = model.components
#     z  = create_z(model,x_c)
#     f! = (F,x) -> Obj_Crit_mix(model, F, 10^x[2], x[1]*prod(model.params.epsilon[i]^z[i] for i in components), x_c)
#     x0 = [1.5,log10(π/6*N_A*sum(z[i]*model.params.segment[i]*model.params.sigma[i]^3 for i in components)/0.15)]
#     r  = nlsolve(f!,x0)
#     T_c = r.zero[1]*prod(model.params.epsilon[i]^z[i] for i in components)
#     v_c = 10^r.zero[2]
#     p_c = get_pressure(model, v_c, T_c, x_c)
#     return (T_c, p_c, v_c)
# end
# #
# function Obj_Crit_mix(model::SAFT, F, v_c,T_c,x_c)
#     fun(x)  = eos(model, create_z(model, [x[1],1-x[1]]), v_c, T_c)
#     df(x)   = ForwardDiff.derivative(fun,x)
#     d2f(x)  = ForwardDiff.derivative(df,x)
#     d3f(x)  = ForwardDiff.derivative(d2f,x)
#     F[1] = d2f(x_c[1])
#     F[2] = d3f(x_c[1])
#     println(F)
#     println(v_c)
#     println(T_c)
# end
## Derivative properties
function get_pressure(model::EoS, v, T, z=[1.]; phase = "unknown")
    z = create_z(model, z)
    fun(x) = eos(model, z, x[1], T)
    df(x)  = ForwardDiff.derivative(fun,x[1])
    return -df(v)
end

function get_entropy(model::EoS, p, T, z=[1.]; phase = "unknown")
    z = create_z(model, z)
    v      = get_volume(model, p, T, z; phase=phase)[1]
    fun(x) = eos(model, z, v, x[1])
    df(x)  = ForwardDiff.derivative(fun,x[1])
    return -df(T)
end

function get_chemical_potential(model::EoS, p, T, z=[1.]; phase = "unknown")
    z = create_z(model, z)
    v      = get_volume(model, p, T, z; phase=phase)[1]
    fun(x) = eos(model, x, v, T)
    df(x)  = ForwardDiff.gradient(fun,x)
    return df(z)
end

function get_internal_energy(model::EoS, p, T, z=[1.]; phase = "unknown")
    z = create_z(model, z)
    v      = get_volume(model, p, T, z; phase=phase)[1]
    fun(x) = eos(model, z, v, x)
    df(x)  = ForwardDiff.derivative(fun,x)
    return fun(T)-df(T)*T
end

function get_enthalpy(model::EoS, p, T, z=[1.]; phase = "unknown")
    z = create_z(model, z)
    v      = get_volume(model, p, T, z; phase=phase)[1]
    fun(x) = eos(model, z, x[1], x[2])
    df(x)  = ForwardDiff.gradient(fun,x)
    return fun([v,T])-df([v,T])[2]*T-df([v,T])[1]*v
end

function get_Gibbs_free_energy(model::EoS, p, T, z=[1.]; phase = "unknown")
    z = create_z(model, z)
    v      = get_volume(model, p, T, z; phase=phase)[1]
    fun(x) = eos(model, z, x[1], T)
    df(x)  = ForwardDiff.derivative(fun,x)
    return fun(v)-df(v)[1]*v
end

function get_Helmholtz_free_energy(model::EoS, p, T, z=[1.]; phase = "unknown")
    z = create_z(model, z)
    v      = get_volume(model, p, T, z; phase=phase)[1]
    fun(x) = eos(model, z, x[1], T)
    df(x)  = ForwardDiff.derivative(fun,x)
    return fun(T)
end

function get_isochoric_heat_capacity(model::EoS, p, T, z=[1.]; phase = "unknown")
    z = create_z(model, z)
    v       = get_volume(model, p, T, z; phase=phase)[1]
    fun(x)  = eos(model, z, v, x)
    df(x)   = ForwardDiff.derivative(fun,x)
    d2f(x)  = ForwardDiff.derivative(df,x)
    return -T*d2f(T)
end

function get_isobaric_heat_capacity(model::EoS, p, T, z=[1.]; phase = "unknown")
    z = create_z(model, z)
    v       = get_volume(model, p, T, z; phase=phase)[1]
    fun(x)  = eos(model, z, x[1], x[2])
    d2f(x)  = ForwardDiff.hessian(fun,x)
    return T*(d2f([v,T])[1,2]^2/d2f([v,T])[1]-d2f([v,T])[2,2])
end

function get_isothermal_compressibility(model::EoS, p, T, z=[1.]; phase = "unknown")
    z = create_z(model, z)
    v       = get_volume(model, p, T, z; phase=phase)[1]
    fun(x)  = eos(model, z, x, T)
    df(x)   = ForwardDiff.derivative(fun,x)
    d2f(x)  = ForwardDiff.derivative(df,x)
    return 1/v*d2f(v)^-1
end

function get_isentropic_compressibility(model::EoS, p, T, z=[1.]; phase = "unknown")
    z = create_z(model, z)
    v       = get_volume(model, p, T, z; phase=phase)[1]
    fun(x)  = eos(model, z, x[1], x[2])
    d2f(x)  = ForwardDiff.hessian(fun,x)
    return 1/v*(d2f([v,T])[1]-d2f([v,T])[1,2]^2/d2f([v,T])[2,2])^-1
end

function get_speed_of_sound(model::EoS, p, T, z=[1.]; phase = "unknown")
    z = create_z(model, z)
    Mr      = sum(z[i]*model.params.Mr[i] for i in model.components)
    v       = get_volume(model, p, T, z; phase=phase)[1]
    fun(x)  = eos(model, z, x[1], x[2])
    d2f(x)  = ForwardDiff.hessian(fun,x)
    return v*sqrt((d2f([v,T])[1]-d2f([v,T])[1,2]^2/d2f([v,T])[2,2])/Mr)
end

function get_isobaric_expansivity(model::EoS, p, T, z=[1.]; phase = "unknown")
    z = create_z(model, z)
    v       = get_volume(model, p, T, z; phase=phase)[1]
    fun(x)  = eos(model, z, x[1], x[2])
    d2f(x)   = ForwardDiff.hessian(fun,x)
    return d2f([v,T])[1,2]/(v*d2f([v,T])[1])
end

function get_Joule_Thomson_coefficient(model::EoS, p, T, z=[1.]; phase = "unknown")
    z = create_z(model, z)
    v       = get_volume(model, p, T, z; phase=phase)[1]
    fun(x)  = eos(model, z, x[1], x[2])
    d2f(x)  = ForwardDiff.hessian(fun,x)
    return -(d2f([v,T])[1,2]-d2f([v,T])[1]*((T*d2f([v,T])[2,2]+v*d2f([v,T])[1,2])/(T*d2f([v,T])[1,2]+v*d2f([v,T])[1])))^-1
end

function get_second_virial_coeff(model::EoS, T, z=[1.])
    V = [1e10]
    z = create_z(model, z)
    fun(x) = eos(model, z, x[1], T)
    df(x)  = ForwardDiff.derivative(fun,x[1])
    d2f(x) = ForwardDiff.derivative(df,x[1])
    return V[1]^2/(R̄*T)*(df(V)+V[1]*d2f(V))
end
