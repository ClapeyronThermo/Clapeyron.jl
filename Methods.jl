using NLopt,NLsolve,DiffResults,ForwardDiff,LinearAlgebra
pai = 3.14159265359
N_A = 6.02214086e23
function Volume(EoS,model,z,p,T,phase="unknown")
    N = length(p)

    ub = [Inf]
    lb = [log10(pai/6*N_A*model.parameters.segment[1]*model.parameters.sigma[1]^3/1)]

    if phase == "unknown" || phase == "liquid"
        x0 = [log10(pai/6*N_A*model.parameters.segment[1]*model.parameters.sigma[1]^3/0.9)]
    elseif phase == "vapour"
        x0 = [log10(pai/6*N_A*model.parameters.segment[1]*model.parameters.sigma[1]^3/1e-2)]
    end

    Vol = []
    if phase == "unknown"
        for i in 1:N
            f(v) = EoS(z[i,:],10^v[1],T[i])+10^v[1]*p[i]
            (f_best,v_best) = Tunneling(f,lb,ub,x0)
            append!(Vol,10^v_best[1])
        end
    else
        opt_min = NLopt.Opt(:LD_MMA, length(ub))
        opt_min.lower_bounds = lb
        opt_min.upper_bounds = ub
        opt_min.xtol_rel = 1e-8
        obj_f0 = x -> f(x)
        obj_f = (x,g) -> NLopt_obj(obj_f0,x,g)
        opt_min.min_objective =  obj_f
        for i in 1:N
            f(v) = EoS(z[i,:],10^v[1],T[i])+10^v[1]*p[i]/T[i]/8.314
            (f_min,v_min) = NLopt.optimize(opt_min, x0)
            append!(Vol,10^v_min[1])
        end
    end
    return Vol
end

function Tunneling(f,lb,ub,x0)
    N = length(ub)
    # Relevant configuration
    tolf=1e-8

    # Minimisation phase
    opt_min = NLopt.Opt(:LD_MMA, length(ub))
    opt_min.lower_bounds = lb
    opt_min.upper_bounds = ub
    opt_min.xtol_rel = 1e-8

    obj_f0 = x -> f(x)
    obj_f = (x,g) -> NLopt_obj(obj_f0,x,g)
    opt_min.min_objective =  obj_f

    (min_f,min_x,status) = NLopt.optimize(opt_min, x0)

    best_f = min_f
    opt_x  = []
    best_x = min_x
    append!(opt_x,[min_x])

    # Tunneling phase
    opt_tun = NLopt.Opt(:LD_MMA, length(ub))
    opt_tun.lower_bounds = lb
    opt_tun.upper_bounds = ub
    opt_tun.xtol_rel = 1e-8
    opt_tun.stopval = -1e-6

    for i in 1:10*N
        T0 = x -> (f(x)-f_best)*prod(exp(1e-2/sqrt(sum((x[i]-x_opt[j][i])^2 for i in 1:N))) for j in 1:k)
        T  = (x,g) -> NLopt_obj(T0,x,g)
        opt_tun.min_objective =  T

        r  = 2.0.*rand(Float64,(N)).-1.0
        ϵ1 = 2*(tolf)^(1/5)*(1+norm(best_x,2))
        x0 = r/norm(r,2).*ϵ1+best_x
        x0 = ub.*(x0.>=ub)+lb.*(x0.<=lb)+x0.*(ub.>x0.>lb)

        # Tunneling
        (new_f,new_x,status) = NLopt.optimize(opt_tun, x0)
        if status != :FORCED_STOP
            println(i)
            break
        end
        # Minimisation
        (min_f,min_x,status) = NLopt.optimize(opt_min, new_x)
        if min_f<best_f
            best_f = min_f
            best_x = min_x
            opt_x  = []
            append!(opt_x,[min_x])
        else
            append!(opt_x,[min_x])
        end

    end
    return (best_f,best_x)
end

function NLopt_obj(f,x,g)
        if length(g) > 0
            df = DiffResults.GradientResult(x)
            df = ForwardDiff.gradient!(df,f,x)
            g .= DiffResults.gradient(df)
            return DiffResults.value(df)
        else
            return f(x)
        end
end

function Psat(EoS,model,T)

    v0    = [log10(pai/6*N_A*model.parameters.segment[1]*model.parameters.sigma[1]^3/0.4),log10(pai/6*N_A*model.parameters.segment[1]*model.parameters.sigma[1]^3/1e-2)]
    v_l   = []
    v_v   = []
    P_sat = []
    for i in 1:length(T)
        f! = (F,x) -> Obj_Sat(F,EoS,model,T[i],10^x[1],10^x[2])
        j! = (J,x) -> Jac_Sat(J,EoS,model,T[i],10^x[1],10^x[2])
        r  =nlsolve(f!,j!,v0)
        append!(v_l,10^r.zero[1])
        append!(v_v,10^r.zero[2])
        append!(P_sat,Pressure(EoS,model,[1],v_v[i],T[i]))
        v0 = r.zero
    end
    return (P_sat,v_l,v_v)
end

function Obj_Sat(F,EoS,model,T,v_l,v_v)
    N_A = 6.02214086e23
    fun(x) = EoS([x[1]],x[2],T)
    df(x)  = ForwardDiff.gradient(fun,x)
    df_l = df([1,v_l[1]])
    df_v = df([1,v_v[1]])
    F[1] = (df_l[2]-df_v[2])*model.parameters.sigma[1]^3*N_A/8.314/model.parameters.epsilon[1]
    F[2] = (df_l[1]-df_v[1])/8.314/model.parameters.epsilon[1]
    println(v_l,v_v)
end

function Jac_Sat(J,EoS,model,T,v_l,v_v)
    N_A = 6.02214086e23
    fun(x) = EoS([x[1]],x[2],T)
    d2f(x) = ForwardDiff.hessian(fun,x)
    d2f_l = d2f([1,v_l[1]])
    d2f_v = d2f([1,v_v[1]])
    J[1] =  v_l[1]*d2f_l[2,2]*model.parameters.sigma[1]^3*N_A*log(10)/8.314/model.parameters.epsilon[1]
    J[1,2] = -v_v[1]*d2f_v[2,2]*model.parameters.sigma[1]^3*N_A*log(10)/8.314/model.parameters.epsilon[1]
    J[2,1] =  v_l[1]*d2f_l[1,2]*log(10)/8.314/model.parameters.epsilon[1]
    J[2,2] = -v_v[1]*d2f_v[1,2]*log(10)/8.314/model.parameters.epsilon[1]
end

function Enthalpy_of_vapourisation(EoS,model,T)
    (P_sat,v_l,v_v) = Psat(EoS,model,T)
    fun(x) = EoS(z,x[1],x[2])
    df(x)  = ForwardDiff.gradient(fun,x)
    H_vap = []
    for i in 1:length(T)
        H_l = fun([v_l[i],T[i]])-df([v_l[i],T[i]])[2]*T[i]-df([v_l[i],T[i]])[1]*v_l[i]
        H_v = fun([v_v[i],T[i]])-df([v_v[i],T[i]])[2]*T[i]-df([v_v[i],T[i]])[1]*v_v[i]
        append!(H_vap,H_v-H_l)
    end
    return H_vap
end

function Pressure(EoS,model,z,p,T,phase="unknown")
    v      = Volume(EoS,model,z,p,T,phase)
    fun(x) = EoS(z,x[1],T)
    df(x)  = ForwardDiff.derivative(fun,x[1])
    return -df(v)
end

function Entropy(EoS,model,z,p,T,phase="unknown")
    v      = Volume(EoS,model,z,p,T,phase)[1]
    fun(x) = EoS(z,v,x[1])
    df(x)  = ForwardDiff.derivative(fun,x[1])
    return -df(T)
end

function Chemical_Potential(EoS,model,z,p,T,phase="unknown")
    v      = Volume(EoS,model,z,p,T,phase)[1]
    fun(x) = EoS(x,v,T)
    df(x)  = ForwardDiff.gradient(fun,x)
    return df(z)
end

function Internal_Energy(EoS,model,z,p,T,phase="unknown")
    v      = Volume(EoS,model,z,p,T,phase)[1]
    fun(x) = EoS(z,v,x)
    df(x)  = ForwardDiff.derivative(fun,x)
    return fun(T)-df(T)*T
end

function Enthalpy(EoS,model,z,p,T,phase="unknown")
    v      = Volume(EoS,model,z,p,T,phase)[1]
    fun(x) = EoS(z,x[1],x[2])
    df(x)  = ForwardDiff.gradient(fun,x)
    return fun([v,T])-df([v,T])[2]*T-df([v,T])[1]*v
end

function Gibbs_free_energy(EoS,model,z,p,T,phase="unknown")
    v      = Volume(EoS,model,z,p,T,phase)[1]
    fun(x) = EoS(z,x[1],T)
    df(x)  = ForwardDiff.derivative(fun,x)
    return fun(v)-df(v)[1]*v
end

function Helmholtz_free_energy(EoS,model,z,p,T,phase="unknown")
    v      = Volume(EoS,model,z,p,T,phase)[1]
    fun(x) = EoS(z,x[1],T)
    df(x)  = ForwardDiff.derivative(fun,x)
    return fun(T)-df(v)[1]*v
end

function Isochoric_heat_capacity(EoS,model,z,p,T,phase="unknown")
    v       = Volume(EoS,model,z,p,T,phase)[1]
    fun(x)  = EoS(z,v,x)
    df(x)   = ForwardDiff.derivative(fun,x)
    d2f(x)  = ForwardDiff.derivative(df,x)
    return -T*d2f(T)
end

function Isobaric_heat_capacity(EoS,model,z,p,T,phase="unknown")
    v       = Volume(EoS,model,z,p,T,phase)[1]
    fun(x)  = EoS(z,x[1],x[2])
    d2f(x)   = ForwardDiff.hessian(fun,x)
    return T*(d2f([v,T])[1,2]^2/d2f([v,T])[1]-d2f([v,T])[2,2])
end

function Isothermal_compressibility(EoS,model,z,p,T,phase="unknown")
    v       = Volume(EoS,model,z,p,T,phase)[1]
    fun(x)  = EoS(z,x,T)
    df(x)   = ForwardDiff.derivative(fun,x)
    d2f(x)  = ForwardDiff.derivative(df,x)
    return 1/v*d2f(v)^-1
end

function Isentropic_compressibility(EoS,model,z,p,T,phase="unknown")
    v       = Volume(EoS,model,z,p,T,phase)[1]
    fun(x)  = EoS(z,x[1],x[2])
    d2f(x)   = ForwardDiff.hessian(fun,x)
    return 1/v*(d2f([v,T])[1]-d2f([v,T])[1,2]^2/d2f([v,T])[2,2])^-1
end

function Speed_of_sound(EoS,model,z,p,T,phase="unknown")
    Mr      = sum(z[i]*model.parameters.Mr[i] for i in model.components)
    v       = Volume(EoS,model,z,p,T,phase)[1]
    fun(x)  = EoS(z,x[1],x[2])
    d2f(x)   = ForwardDiff.hessian(fun,x)
    return v*sqrt((d2f([v,T])[1]-d2f([v,T])[1,2]^2/d2f([v,T])[2,2])/Mr)
end

function Isobaric_expansivity(EoS,model,z,p,T,phase="unknown")
    v       = Volume(EoS,model,z,p,T,phase)[1]
    fun(x)  = EoS(z,x[1],x[2])
    d2f(x)   = ForwardDiff.hessian(fun,x)
    return d2f([v,T])[1,2]/(v*d2f([v,T])[1])
end

function Joule_Thomson_coefficient(EoS,model,z,p,T,phase="unknown")
    v       = Volume(EoS,model,z,p,T,phase)[1]
    fun(x)  = EoS(z,x[1],x[2])
    d2f(x)   = ForwardDiff.hessian(fun,x)
    return -(d2f([v,T])[1,2]-d2f([v,T])[1]*((T*d2f([v,T])[2,2]+v*d2f([v,T])[1,2])/(T*d2f([v,T])[1,2]+v*d2f([v,T])[1])))^-1
end
