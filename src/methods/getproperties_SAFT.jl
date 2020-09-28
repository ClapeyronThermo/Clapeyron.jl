## Standard pressure solver
function get_volume(model::SAFT, z, p, T, phase="unknown")
    components = model.components

    N = length(p)

    ub = [Inf]
    lb = [log10(π/6*N_A*sum(z[i]*model.params.segment[i]*model.params.sigma[i]^3 for i in model.components)/1)]

    if phase == "unknown" || phase == "liquid"
        x0 = [log10(π/6*N_A*sum(z[i]*model.params.segment[i]*model.params.sigma[i]^3 for i in model.components)/0.5)]
    elseif phase == "vapour"
        x0 = [log10(π/6*N_A*sum(z[i]*model.params.segment[i]*model.params.sigma[i]^3 for i in model.components)/1e-2)]
    end

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


## Saturation conditions solver
function get_Psat(model::SAFT, T)
    components = model.components
    v0    = [log10(π/6*N_A*model.params.segment[components[1]]*model.params.sigma[components[1]]^3/0.4),
             log10(π/6*N_A*model.params.segment[components[1]]*model.params.sigma[components[1]]^3/1e-3)]
    v_l   = []
    v_v   = []
    P_sat = []
    for i in 1:length(T)
        f! = (F,x) -> Obj_Sat(model, F, T[i], 10^x[1], 10^x[2])
        j! = (J,x) -> Jac_Sat(model, J, T[i], 10^x[1], 10^x[2])
        r  =nlsolve(f!,j!,v0)
        append!(v_l,10^r.zero[1])
        append!(v_v,10^r.zero[2])
        append!(P_sat,get_pressure(model,create_z(model, [1.0]),v_v[i],T[i]))
        v0 = r.zero
    end
    return (P_sat,v_l,v_v)
end

function Obj_Sat(model::SAFT, F, T, v_l, v_v)
    components = model.components
    fun(x) = eos(model, create_z(model, [x[1]]), x[2], T)
    df(x)  = ForwardDiff.gradient(fun,x)
    df_l = df([1,v_l[1]])
    df_v = df([1,v_v[1]])
    F[1] = (df_l[2]-df_v[2])*model.params.sigma[components[1]]^3*N_A/R̄/model.params.epsilon[components[1]]
    F[2] = (df_l[1]-df_v[1])/R̄/model.params.epsilon[components[1]]
end

function Jac_Sat(model::SAFT, J, T, v_l, v_v)
    components = model.components
    fun(x) = eos(model, create_z(model, [x[1]]), x[2], T)
    d2f(x) = ForwardDiff.hessian(fun,x)
    d2f_l = d2f([1,v_l[1]])
    d2f_v = d2f([1,v_v[1]])
    J[1] =  v_l[1]*d2f_l[2,2]*model.params.sigma[components[1]]^3*N_A*log(10)/R̄/model.params.epsilon[components[1]]
    J[1,2] = -v_v[1]*d2f_v[2,2]*model.params.sigma[components[1]]^3*N_A*log(10)/R̄/model.params.epsilon[components[1]]
    J[2,1] =  v_l[1]*d2f_l[1,2]*log(10)/R̄/model.params.epsilon[components[1]]
    J[2,2] = -v_v[1]*d2f_v[1,2]*log(10)/R̄/model.params.epsilon[components[1]]
end

## Critical point solver
function get_Pcrit(model::SAFT)
    components = model.components
    f! = (F,x) -> Obj_Crit(model, F, x[1]*model.params.epsilon[components[1]], 10^x[2])
    # j! = (J,x) -> Jac_Crit(J,eos,model,x[1]*model.params.epsilon[(1, 1)],10^x[2])
    x0 = [1.5, log10(π/6*N_A*model.params.segment[components[1]]*model.params.sigma[components[1]]^3/0.3)]
    r  = nlsolve(f!,x0)
    T_c = r.zero[1]*model.params.epsilon[components[1]]
    v_c = 10^r.zero[2]
    p_c = get_pressure(model, create_z(model, [1.0]), v_c, T_c)
    return (T_c, p_c, v_c)
end

function Obj_Crit(model::SAFT, F, T_c, v_c)
    fun(x)  = eos(model, create_z(model, [1]), x[1], T_c)
    df(x)   = ForwardDiff.derivative(fun,x)
    d2f(x)  = ForwardDiff.derivative(df,x)
    d3f(x)  = ForwardDiff.derivative(d2f,x)
    F[1] = d2f(v_c)
    F[2] = d3f(v_c)
end

# function Jac_Crit(F,eos,model,T_c,v_c)
#     fun(x)  = eos([1],x[1],x[2])
#     df(x)   = ForwardDiff.gradient(fun,x)
#     d2f(x)  = ForwardDiff.gradient(df,x)
#     d3f(x)  = ForwardDiff.gradient(d2f,x)
#     d4f(x)  = ForwardDiff.gradient(d3f,x)
#
#     F[1] = d2f(v_c)
#     F[2] = d3f(v_c)
# end

function get_enthalpy_vap(model::SAFT, T)
    (P_sat,v_l,v_v) = get_Psat(model,T)
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

## Derivative properties
function get_pressure(model::SAFT, z, v, T, phase="unknown")
    fun(x) = eos(model, z, x[1], T)
    df(x)  = ForwardDiff.derivative(fun,x[1])
    return -df(v)
end

function get_entropy(model::SAFT, z, p, T, phase="unknown")
    v      = get_volume(model, z, p, T, phase)[1]
    fun(x) = eos(model, z, v, x[1])
    df(x)  = ForwardDiff.derivative(fun,x[1])
    return -df(T)
end

function get_chemical_potential(model::SAFT, z, p, T, phase="unknown")
    v      = get_volume(model, z, p, T, phase)[1]
    fun(x) = eos(model, x, v, T)
    df(x)  = ForwardDiff.gradient(fun,x)
    return df(z)
end

function get_internal_energy(model::SAFT, z, p, T, phase="unknown")
    v      = get_volume(model, z, p, T, phase)[1]
    fun(x) = eos(model, z, v, x)
    df(x)  = ForwardDiff.derivative(fun,x)
    return fun(T)-df(T)*T
end

function get_enthalpy(model::SAFT, z, p, T, phase="unknown")
    v      = get_volume(model, z, p, T, phase)[1]
    fun(x) = eos(model, z, x[1], x[2])
    df(x)  = ForwardDiff.gradient(fun,x)
    return fun([v,T])-df([v,T])[2]*T-df([v,T])[1]*v
end

function get_Gibbs_free_energy(model::SAFT, z, p, T, phase="unknown")
    v      = get_volume(model, z, p, T, phase)[1]
    fun(x) = eos(model, z, x[1], T)
    df(x)  = ForwardDiff.derivative(fun,x)
    return fun(v)-df(v)[1]*v
end

function get_Helmholtz_free_energy(model::SAFT, z, p, T, phase="unknown")
    v      = get_volume(model, z, p, T, phase)[1]
    fun(x) = eos(model, z, x[1], T)
    df(x)  = ForwardDiff.derivative(fun,x)
    return fun(T)-df(v)[1]*v
end

function get_isochoric_heat_capacity(model::SAFT, z, p, T, phase="unknown")
    v       = get_volume(model, z, p, T, phase)[1]
    fun(x)  = eos(model, z, v, x)
    df(x)   = ForwardDiff.derivative(fun,x)
    d2f(x)  = ForwardDiff.derivative(df,x)
    return -T*d2f(T)
end

function get_isobaric_heat_capacity(model::SAFT, z, p, T, phase="unknown")
    v       = get_volume(model, z, p, T, phase)[1]
    fun(x)  = eos(model, z, x[1], x[2])
    d2f(x)  = ForwardDiff.hessian(fun,x)
    return T*(d2f([v,T])[1,2]^2/d2f([v,T])[1]-d2f([v,T])[2,2])
end

function get_isothermal_compressibility(model::SAFT, z, p, T, phase="unknown")
    v       = get_volume(model, z, p, T, phase)[1]
    fun(x)  = eos(model, z, x, T)
    df(x)   = ForwardDiff.derivative(fun,x)
    d2f(x)  = ForwardDiff.derivative(df,x)
    return 1/v*d2f(v)^-1
end

function get_isentropic_compressibility(model::SAFT, z, p, T, phase="unknown")
    v       = get_volume(model, z, p, T, phase)[1]
    fun(x)  = eos(model, z, x[1], x[2])
    d2f(x)  = ForwardDiff.hessian(fun,x)
    return 1/v*(d2f([v,T])[1]-d2f([v,T])[1,2]^2/d2f([v,T])[2,2])^-1
end

function get_speed_of_sound(model::SAFT, z, p, T, phase="unknown")
    Mr      = sum(z[i]*model.params.Mr[i] for i in model.components)
    v       = get_volume(model, z, p, T, phase)[1]
    fun(x)  = eos(model, z, x[1], x[2])
    d2f(x)  = ForwardDiff.hessian(fun,x)
    return v*sqrt((d2f([v,T])[1]-d2f([v,T])[1,2]^2/d2f([v,T])[2,2])/Mr)
end

function get_isobaric_expansivity(model::SAFT, z, p, T, phase="unknown")
    v       = get_volume(model, z, p, T, phase)[1]
    fun(x)  = eos(model, z, x[1], x[2])
    d2f(x)   = ForwardDiff.hessian(fun,x)
    return d2f([v,T])[1,2]/(v*d2f([v,T])[1])
end

function get_Joule_Thomson_coefficient(model::SAFT, z, p, T, phase="unknown")
    v       = get_volume(model, z, p, T, phase)[1]
    fun(x)  = eos(model, z, x[1], x[2])
    d2f(x)  = ForwardDiff.hessian(fun,x)
    return -(d2f([v,T])[1,2]-d2f([v,T])[1]*((T*d2f([v,T])[2,2]+v*d2f([v,T])[1,2])/(T*d2f([v,T])[1,2]+v*d2f([v,T])[1])))^-1
end
