#derivative logic

#ForwardDiff compiles one method per function, so separating the
#differentiation logic from the property logic allows the differentials
#to be compiled only once

function ∂f∂t(model,v,t,z)
    return ForwardDiff.derivative(∂t -> eos(model,z,v,∂t),t)
end

function ∂f∂v(model,v,t,z)
    return ForwardDiff.derivative(∂v -> eos(model,z,∂v,t),v)
end

#returns a tuple of the form ([∂f∂v,∂f∂t],f),using the least amount of computation
function ∂f(model,v,t,z)
    f(w) = eos(model,z,first(w),last(w))
    v,t = promote(v,t)
    vt_vec = SVector(v,t)
    ∂result = DiffResults.GradientResult(vt_vec)
    res_∂f =  ForwardDiff.gradient!(∂result, f,vt_vec)
    _f =  DiffResults.value(res_∂f)
    _∂f = DiffResults.gradient(res_∂f)
    return (_∂f,_f)
end

#returns a tuple, of the form (hess_vt(f),grad_vt(f),f), it does one allocation because of a bug
function ∂2f(model,v,t,z)
    f(w) = eos(model,z,first(w),last(w))
    v,t = promote(v,t)
    vt_vec =   SVector(v,t)
    ∂result = DiffResults.HessianResult(vt_vec)
    res_∂f =  ForwardDiff.hessian!(∂result, f,vt_vec)
    _f =  DiffResults.value(res_∂f)
    _∂f = DiffResults.gradient(res_∂f)
    _∂2f = DiffResults.hessian(res_∂f)
    return (_∂2f,_∂f,_f)
end

#returns hess_vt of helmholtz energy
function f_hess(model,v,t,z)
    f(w) = eos(model,z,first(w),last(w))
    v,t = promote(v,t)
    vt_vec = SVector(v,t)
    return ForwardDiff.hessian(f,vt_vec)
end

## Standard pressure solver

function get_volume(model::EoS, p, T, z=[1.]; phase = "unknown")
    z = create_z(model, z)
    components = model.components

    N = length(p)

    ub = [Inf]
    lb = lb_volume(model,z; phase = phase)
    x0 = x0_volume(model,z; phase = phase)
    f = v -> eos(model, z, exp10(v[1]), T) + exp10(v[1])*p
    if phase == "unknown"
        (f_best,v_best) = Solvers.tunneling(f,lb,ub,x0)
        return exp10(v_best[1])
    else
        opt_min = NLopt.Opt(:LD_MMA, length(ub))
        opt_min.lower_bounds = lb
        opt_min.upper_bounds = ub
        opt_min.xtol_rel     = 1e-8
        obj_f0 = x -> f(x)
        obj_f  = (x,g) -> Solvers.NLopt_obj(obj_f0,x,g)
        opt_min.min_objective =  obj_f
        (f_min,v_min) = NLopt.optimize(opt_min, x0)
        return exp10(v_min[1])
    end
end


## Pure saturation conditions solver
function get_sat_pure(model::EoS, T; v0 = nothing)
    components = model.components

    f! = (F,x) -> Obj_Sat(model, F, T, exp10(x[1]), exp10(x[2]))
    j! = (J,x) -> Jac_Sat(model, J, T, exp10(x[1]), exp10(x[2]))
    fj! = (F,J,x) -> Obj_Jac_sat(model,F,J,T,exp10(x[1]), exp10(x[2]))
    jv! = (x) -> Jvop_sat(x, model, T)
    vectorobj = NLSolvers.VectorObjective(f!,j!,fj!,jv!)
    vectorprob = NEqProblem(vectorobj)
    if v0 == nothing
        try
            v0    = x0_sat_pure(model)
            (P_sat,v_l,v_v) = solve_sat_pure(model,v0,vectorprob,T)
            if abs(v_l-v_v)/v_l<1e-2
                v0    = x0_sat_pure(model)
                v0[1] = v0[1]+log10(0.5/0.52)
                (P_sat,v_l,v_v) = solve_sat_pure(model,v0,vectorprob,T)
                if abs(v_l-v_v)/v_l<1e-2
                    v0    = x0_sat_pure(model)
                    v0[1] = v0[1]+log10(0.5/0.48)
                    (P_sat,v_l,v_v) = solve_sat_pure(model,v0,vectorprob,T)
                end
                return (P_sat,v_l,v_v)

            end
            return (P_sat,v_l,v_v)
        catch y
            if isa(y, DomainError)
                try
                    v0    = x0_sat_pure(model)
                    v0[1] = v0[1]+log10(0.5/0.3)
                    (P_sat,v_l,v_v) = solve_sat_pure(model,v0,vectorprob,T)
                    if abs(v_l-v_v)/v_l<1e-2
                        v0    = x0_sat_pure(model)
                        v0[1] = v0[1]+log10(0.5/0.32)
                        (P_sat,v_l,v_v) = solve_sat_pure(model,v0,vectorprob,T)
                        if abs(v_l-v_v)/v_l<1e-2
                            v0    = x0_sat_pure(model)
                            v0[1] = v0[1]+log10(0.5/0.28)
                            (P_sat,v_l,v_v) = solve_sat_pure(model,v0,vectorprob,T)
                        end
                    end
                    return (P_sat,v_l,v_v)
                catch y
                    v0    = x0_sat_pure(model)
                    v0[1] = v0[1]+log10(0.5/0.4)
                    (P_sat,v_l,v_v) = solve_sat_pure(model,v0,vectorprob,T)
                    if abs(v_l-v_v)/v_l<1e-2
                        v0    = x0_sat_pure(model)
                        v0[1] = v0[1]+log10(0.5/0.42)
                        (P_sat,v_l,v_v) = solve_sat_pure(model,v0,vectorprob,T)
                        if abs(v_l-v_v)/v_l<1e-2
                            v0    = x0_sat_pure(model)
                            v0[1] = v0[1]+log10(0.5/0.38)
                            (P_sat,v_l,v_v) = solve_sat_pure(model,v0,vectorprob,T)
                        end
                    end
                    return (P_sat,v_l,v_v)
                end

            end
        end
    else
        (P_sat,v_l,v_v) = solve_sat_pure(model,v0,vectorprob,T)
        return (P_sat,v_l,v_v)
    end

end

function solve_sat_pure(model::EoS,v0,vectorprob,T)
    r = solve(vectorprob, v0, LineSearch(Newton()), NEqOptions())
    v_l = exp10(r.info.solution[1])
    v_v = exp10(r.info.solution[2])
    P_sat = get_pressure(model,exp10(r.info.solution[2]),T)
    return (P_sat,v_l,v_v)
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

function Obj_Jac_sat(model::EoS, F, J, T, v_l, v_v)
    Obj_Sat(model, F, T, v_l, v_v)
    Jac_Sat(model, J, T, v_l, v_v)
    F, J
end

function Jvop_sat(x,model::EoS,T)
    function Jac_satV(Fv, v)
        components = model.components
        fun(z) = eos(model, create_z(model, [z[1]]), z[2], T)
        d2f(z) = ForwardDiff.hessian(fun,z)
        d2f_l = d2f([1,x[1]])
        d2f_v = d2f([1,x[2]])
        (p_scale,μ_scale) = scale_sat_pure(model)
        Fv[1,] =  x[1]*d2f_l[2,2]*log(10)*p_scale*v[1]-x[2]*d2f_v[2,2]*log(10)*p_scale*v[2]
        Fv[2,] =  x[1]*d2f_l[1,2]*log(10)*μ_scale*v[1]-x[2]*d2f_v[1,2]*log(10)*μ_scale*v[2]
    end
    LinearMap(Jac_satV, length(x))
end

function get_enthalpy_vap(model::EoS, T)
    (P_sat,v_l,v_v) = get_sat_pure(model,T)
    fun(x) = eos(model, create_z(model,[1.0]), x[1], x[2])
    df(x)  = ForwardDiff.gradient(fun,x)
    H_l = fun([v_l,T])-df([v_l,T])[2]*T-df([v_l,T])[1]*v_l
    H_v = fun([v_v,T])-df([v_v,T])[2]*T-df([v_v,T])[1]*v_v
    H_vap=H_v-H_l
    return H_vap
end
## Pure critical point solver
function get_crit_pure(model::EoS; units = false, output=[u"K", u"Pa", u"m^3"])
    components = model.components
    T̄  = T_crit_pure(model)
    f! = (F,x) -> Obj_Crit(model, F, x[1]*T̄, exp10(x[2]))
    # j! = (J,x) -> Jac_Crit(J,eos,model,x[1]*model.params.epsilon[(1, 1)],exp10(x[2]))
    x0 = x0_crit_pure(model)
    r  = nlsolve(f!,x0)
    T_c = r.zero[1]*T̄
    v_c = exp10(r.zero[2])
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
function get_bubble_pressure(model, T, x; v0 =nothing)
    components = model.components
    if v0 == nothing
        y0    = 10 .*x[1,:]./(1 .+x[1,:].*(10-1))
        y0    = y0 ./sum(y0[i] for i in 1:length(components))
        X     = create_z(model,x[1,:])
        Y0    = create_z(model,y0)
        v0    = [log10(π/6*N_A*sum(X[i]*model.params.segment[i]*model.params.sigma[i]^3 for i in components)/0.45),
                 log10(π/6*N_A*sum(Y0[i]*model.params.segment[i]*model.params.sigma[i]^3 for i in components)/1e-4)]
        append!(v0,y0[1:end-1])
    end
    v_l   = []
    v_v   = []
    y     = deepcopy(x)
    P_sat = []
    for i in 1:size(x)[1]
        f! = (F,z) -> Obj_bubble_pressure(model, F, T, exp10(z[1]), exp10(z[2]), x[i,:], z[3:end])
        j! = (J,z) -> Jac_bubble_pressure(model, J, T, exp10(z[1]), exp10(z[2]), x[i,:], z[3:end])
        r  =nlsolve(f!,j!,v0)
        append!(v_l,exp10(r.zero[1]))
        append!(v_v,exp10(r.zero[2]))
        append!(P_sat,get_pressure(model,v_l[i],T,x[i,:]))
        y[i,1:end-1] = r.zero[3:end]
        y[i,end] = 1-sum(r.zero[3:end])
        v0 = r.zero
    end
    return (P_sat, v_l, v_v, y)
end

function Obj_bubble_pressure(model, F, T, v_l, v_v, x, y)
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

function Jac_bubble_pressure(model, J, T, v_l, v_v, x, y)
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
#     f! = (F,x) -> Obj_Crit_mix(model, F, exp10(x[2]), x[1]*prod(model.params.epsilon[i]^z[i] for i in components), x_c)
#     x0 = [1.5,log10(π/6*N_A*sum(z[i]*model.params.segment[i]*model.params.sigma[i]^3 for i in components)/0.15)]
#     r  = nlsolve(f!,x0)
#     T_c = r.zero[1]*prod(model.params.epsilon[i]^z[i] for i in components)
#     v_c = exp10(r.zero[2])
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
    return -∂f∂v(model,v,T,z)
end

function get_entropy(model::EoS, p, T, z=[1.]; phase = "unknown")
    z = create_z(model, z)
    v      = get_volume(model, p, T, z; phase=phase)[1]
    return -∂f∂t(model,v,T,z)
end

function get_chemical_potential(model::EoS, p, T, z=[1.]; phase = "unknown")
    z = create_z(model, z)
    v      = get_volume(model, p, T, z; phase=phase)
    fun(x) = eos(model, x, v, T)
    return ForwardDiff.gradient(fun,z)
end

function get_internal_energy(model::EoS, p, T, z=[1.]; phase = "unknown")
    z = create_z(model, z)
    v      = get_volume(model, p, T, z; phase=phase)[1]
    _df,f =  ∂f(model,v,T,z)
    dv,dt = _df
    return f  - dt*T
end

function get_enthalpy(model::EoS, p, T, z=[1.]; phase = "unknown")
    z = create_z(model, z)
    v      = get_volume(model, p, T, z; phase=phase)
    _df,f =  ∂f(model,v,T,z)
    dv,dt = _df
    return f  - dv*v - dt*T
end

function get_Gibbs_free_energy(model::EoS, p, T, z=[1.]; phase = "unknown")
    z = create_z(model, z)
    v      = get_volume(model, p, T, z; phase=phase)
    _df,f =  ∂f(model,v,T,z)
    dv,dt = _df
    return f  - dv*v
end

function get_Helmholtz_free_energy(model::EoS, p, T, z=[1.]; phase = "unknown")
    z = create_z(model, z)
    v      = get_volume(model, p, T, z; phase=phase)
    return eos(model, z, v, T)

end

function get_isochoric_heat_capacity(model::EoS, p, T, z=[1.]; phase = "unknown")
    z = create_z(model, z)
    v       = get_volume(model, p, T, z; phase=phase)
    fun(x)  = eos(model, z, v, x)
    df(x)   = ForwardDiff.derivative(fun,x)
    d2f(x)  = ForwardDiff.derivative(df,x)
    return -T*d2f(T)
end

function get_isobaric_heat_capacity(model::EoS, p, T, z=[1.]; phase = "unknown")
    z = create_z(model, z)
    v       = get_volume(model, p, T, z; phase=phase)
    d2f = f_hess(model,v,T,z)
    return T*(d2f[1,2]^2/d2f[1]-d2f[2,2])
end

function get_isothermal_compressibility(model::EoS, p, T, z=[1.]; phase = "unknown")
    z = create_z(model, z)
    v       = get_volume(model, p, T, z; phase=phase)
    d2f = f_hess(model,v,T,z)
    return 1/v*d2f[1,1]^-1
end

function get_isentropic_compressibility(model::EoS, p, T, z=[1.]; phase = "unknown")
    z = create_z(model, z)
    v       = get_volume(model, p, T, z; phase=phase)
    d2f = f_hess(model,v,T,z)
    return 1/v*(d2f[1]-d2f[1,2]^2/d2f[2,2])^-1
end

function get_speed_of_sound(model::EoS, p, T, z=[1.]; phase = "unknown")
    z = create_z(model, z)
    Mr      = sum(z[i]*model.params.Mr[i] for i in model.components)
    v       = get_volume(model, p, T, z; phase=phase)[1]
    d2f = f_hess(model,v,T,z)
    return v*sqrt((d2f[1]-d2f[1,2]^2/d2f[2,2])/Mr)
end

function get_isobaric_expansivity(model::EoS, p, T, z=[1.]; phase = "unknown")
    z = create_z(model, z)
    v = get_volume(model, p, T, z; phase=phase)
    d2f = f_hess(model,v,T,z)
    return d2f[1,2]/(v*d2f[1])
end

function get_Joule_Thomson_coefficient(model::EoS, p, T, z=[1.]; phase = "unknown")
    z = create_z(model, z)
    v  = get_volume(model, p, T, z; phase=phase)
    d2f = f_hess(model,v,T,z)
    return -(d2f[1,2]-d2f[1]*((T*d2f[2,2]+v*d2f[1,2])/(T*d2f[1,2]+v*d2f[1])))^-1
end

function get_second_virial_coeff(model::EoS, T, z=[1.])
    V = 1e10
    z = create_z(model, z)
    _∂2f = ∂2f(model,V,T,z)
    hessf,gradf,f = _∂2f
    return V^2/(R̄*T)*(gradf[1]+V*hessf[1])
end
