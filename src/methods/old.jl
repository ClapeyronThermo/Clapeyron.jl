## Old pressure solver


function volume(model::EoSModel, p, T,  z=SA[1.]; phase = "unknown")
    N = length(p)

    ub = [Inf]
    lb = lb_volume(model,z; phase = phase)
    
    x0 = x0_volume(model,p,T,z; phase = phase)
    f = v -> eos(model, exp10(v[1]), T,z) + exp10(v[1])*p
    #looking at best phase using tunneling
    if phase == "unknown"
        (f_best,v_best) = Solvers.tunneling(f,lb,ub,x0)
        #@show eos.(Ref(model),exp10.(v_best),T,Ref(z)) .- p.*exp10.(v_best)

        return exp10(v_best[1])
    else
        _ub = 100*T/p
        _lb = exp10(only(lb))
        f0(vx) = pressure(model,vx,T,z,phase=phase) - p        
        _v0 = exp10(only(x0))
        #try direct newton solving
        vobj = Solvers.ad_newton(f0,_v0)
        if (_lb <= vobj <= _ub)
            return vobj
        else
            vobj = Roots.find_zero(f0,(_lb,_ub),FalsePosition())
        end
        #=
        opt_min = NLopt.Opt(:LD_MMA, length(ub))
        opt_min.lower_bounds = lb
        opt_min.upper_bounds = ub
        opt_min.xtol_rel     = 1e-8
        obj_f0 = x -> f(x)
        obj_f  = (x,g) -> Solvers.NLopt_obj(obj_f0,x,g)
        opt_min.min_objective =  obj_f
        (f_min,v_min) = NLopt.optimize(opt_min, x0)
        #@show eos.(Ref(model),exp10.(v_min),T,Ref(z)) .- p.*exp10.(v_min)
        return exp10(v_min[1])
        =#
    end
end

function Jac_Sat(model::EoSModel, J, T, v_l, v_v)
    #components = model.components
    fun(x) = eos(model, x[2], T,SA[x[1]])
    d2f(x) = ForwardDiff.hessian(fun,x)
    d2f_l = d2f(SA[one(v_l),v_l])
    d2f_v = d2f(SA[one(v_v),v_v])
    (p_scale,μ_scale) = scale_sat_pure(model)
    J[1,1] =  v_l[1]*d2f_l[2,2]*log(10)*p_scale
    J[1,2] = -v_v[1]*d2f_v[2,2]*log(10)*p_scale
    J[2,1] =  v_l[1]*d2f_l[1,2]*log(10)*μ_scale
    J[2,2] = -v_v[1]*d2f_v[1,2]*log(10)*μ_scale
end

function Obj_Jac_sat(model::EoSModel, F, J, T, v_l, v_v)
    Obj_Sat(model, F, T, v_l, v_v)
    Jac_Sat(model, J, T, v_l, v_v)
    F, J
end

function Jvop_sat(x,model::EoSModel,T)
    function Jac_satV(Fv, v)
        #components = model.components
        fun(x) = eos(model, x[2], T,[x[1]])
        d2f(z) = ForwardDiff.hessian(fun,z)
        d2f_l = d2f([1,x[1]])
        d2f_v = d2f([1,x[2]])
        (p_scale,μ_scale) = scale_sat_pure(model)
        Fv[1,] =  x[1]*d2f_l[2,2]*log(10)*p_scale*v[1]-x[2]*d2f_v[2,2]*log(10)*p_scale*v[2]
        Fv[2,] =  x[1]*d2f_l[1,2]*log(10)*μ_scale*v[1]-x[2]*d2f_v[1,2]*log(10)*μ_scale*v[2]
    end
    LinearMap(Jac_satV, length(x))
end




## Mixture critical point solver
# function crit_mix(model::SAFT,x_c)
#     components = model.components
#     z  = create_z(model,x_c)
#     f! = (F,x) -> Obj_Crit_mix(model, F, exp10(x[2]), x[1]*prod(model.params.epsilon[i]^z[i] for i in components), x_c)
#     x0 = [1.5,log10(π/6*N_A*sum(z[i]*model.params.segment[i]*model.params.sigma[i]^3 for i in components)/0.15)]
#     r  = Solvers.nlsolve(f!,x0)
#     T_c = r.zero[1]*prod(model.params.epsilon[i]^z[i] for i in components)
#     v_c = exp10(r.zero[2])
#     p_c = pressure(model, v_c, T_c, x_c)
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