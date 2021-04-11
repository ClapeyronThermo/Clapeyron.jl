
function volume(model::EoSModel,p,T,z=SA[1.0];phase="unknown")

    fp(_v) = pressure(model,_v,T,z) - p

#Threaded version

    if is_liquid(phase)
        return volume_compress(model,p,T,z)
    elseif is_vapour(phase) | is_supercritical(phase)
        vg0 = volume_virial(model,p,T,z)
        return Solvers.ad_newton(fp,vg0,rtol=1e-08)
    else #unknown handled here
    end
    
    _vg = Threads.@spawn begin
        vg0 = volume_virial(model,$p,$T,$z)
        Solvers.ad_newton(fp,vg0,rtol=1e-08)
    end
        
    _vl = Threads.@spawn volume_compress(model,$p,$T,$z)
    #fp(_v) = pressure(model,_v,T,z) - p
    vg = fetch(_vg)
    vl = fetch(_vl)


# Serial version
#=
    vg0 = volume_virial(model,p,T,z)
    vg =Solvers.ad_newton(fp,vg0,rtol=1e-08)

    vl =  volume_compress(model,p,T,z)
=#
    function gibbs(v)
        _df,f =  ∂f(model,v,T,z)
        dv,dt = _df
        if abs((p+dv)/p) > 0.03
            return Inf
        else
            return f  +p*v
        end
    end
    #this catches the supercritical phase as well
    if vl ≈ vg
        return vl
    end  
        gg = gibbs(vg)
        gl = gibbs(vl)
        #@show vg,vl
        #@show gg,gl
        return ifelse(gg<gl,vg,vl)
    
end

## Old pressure solver

#=
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
=#


## Pure saturation conditions solver
function sat_pure(model::EoSModel, T; v0 = nothing)
    #components = model.components
    f! = (F,x) -> Obj_Sat(model, F, T, exp10(x[1]), exp10(x[2]))
    #j! = (J,x) -> Jac_Sat(model, J, T, exp10(x[1]), exp10(x[2]))
    #fj! = (F,J,x) -> Obj_Jac_sat(model,F,J,T,exp10(x[1]), exp10(x[2]))
    #jv! = (x) -> Jvop_sat(x, model, T)
    #vectorobj = NLSolvers.VectorObjective(f!,j!,fj!,jv!)
    #vectorprob = NEqProblem(vectorobj)
    #vectorprob = f!
    if v0 === nothing
        try
            v0    = x0_sat_pure(model,T)
            (P_sat,v_l,v_v) = sat_pure(model,v0,f!,T)
            if abs(v_l-v_v)/v_l<1e-2
                v0    = x0_sat_pure(model,T)
                v0[1] = v0[1]+log10(0.5/0.52)
                (P_sat,v_l,v_v) = sat_pure(model,v0,f!,T)
                if abs(v_l-v_v)/v_l<1e-2
                    v0    = x0_sat_pure(model,T)
                    v0[1] = v0[1]+log10(0.5/0.48)
                    (P_sat,v_l,v_v) = sat_pure(model,v0,f!,T)
                end
                return (P_sat,v_l,v_v)

            end
            return (P_sat,v_l,v_v)
        catch y
            if isa(y, DomainError)
                try
                    v0    = x0_sat_pure(model,T)
                    v0[1] = v0[1]+log10(0.5/0.3)
                    (P_sat,v_l,v_v) = sat_pure(model,v0,f!,T)
                    if abs(v_l-v_v)/v_l<1e-2
                        v0    = x0_sat_pure(model,T)
                        v0[1] = v0[1]+log10(0.5/0.32)
                        (P_sat,v_l,v_v) = sat_pure(model,v0,f!,T)
                        if abs(v_l-v_v)/v_l<1e-2
                            v0    = x0_sat_pure(model,T)
                            v0[1] = v0[1]+log10(0.5/0.28)
                            (P_sat,v_l,v_v) = sat_pure(model,v0,f!,T)
                        end
                    end
                    return (P_sat,v_l,v_v)
                catch y
                    v0    = x0_sat_pure(model,T)
                    v0[1] = v0[1]+log10(0.5/0.4)
                    (P_sat,v_l,v_v) = sat_pure(model,v0,f!,T)
                    if abs(v_l-v_v)/v_l<1e-2
                        v0    = x0_sat_pure(model,T)
                        v0[1] = v0[1]+log10(0.5/0.42)
                        (P_sat,v_l,v_v) = sat_pure(model,v0,f!,T)
                        if abs(v_l-v_v)/v_l<1e-2
                            v0    = x0_sat_pure(model,T)
                            v0[1] = v0[1]+log10(0.5/0.38)
                            (P_sat,v_l,v_v) = sat_pure(model,v0,f!,T)
                        end
                    end
                    return (P_sat,v_l,v_v)
                end
            else
                rethrow(y)
            end
        end
    else
        (P_sat,v_l,v_v) = sat_pure(model,v0,vectorprob,T)
        return (P_sat,v_l,v_v)
    end
end

function sat_pure(model::EoSModel,v0,f!,T)
    #@show v0
    r = Solvers.nlsolve(f!, v0)
    v_l = exp10(r.info.solution[1])
    v_v = exp10(r.info.solution[2])
    P_sat = pressure(model,exp10(r.info.solution[2]),T)
    return (P_sat,v_l,v_v)
end

function Obj_Sat(model::EoSModel, F, T, v_l, v_v)
    #components = model.components
    fun(x) = eos(model, x[2], T,SA[x[1]])
    df(x)  = ForwardDiff.gradient(fun,x)
    df_l = df(SA[one(v_l),v_l])
    df_v = df(SA[one(v_v),v_v])
    (p_scale,μ_scale) = scale_sat_pure(model)
    F[1] = (df_l[2]-df_v[2])*p_scale
    F[2] = (df_l[1]-df_v[1])*μ_scale
end
#=
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
=#

function enthalpy_vap(model::EoSModel, T)
    (P_sat,v_l,v_v) = sat_pure(model,T)
   #= _dfl,fl =  ∂f(model,v_l,T,SA[1.0])
    _dfv,fv =  ∂f(model,v_v,T,SA[1.0])
    dvl,dtl = _dfl
    dvv,dtv = _dfv
    H_l = fl  - dvl*v_l - dtl*T
    H_v = fv  - dvv*v_v - dtv*T =#
    H_v = vt_enthalpy(model,v_v,T)
    H_l = vt_enthalpy(model,v_l,T)
    H_vap=H_v-H_l
    return H_vap
end
## Pure critical point solver
function crit_pure(model::EoSModel)
    T̄  = T_crit_pure(model)
    f! = (F,x) -> Obj_Crit(model, F, x[1]*T̄, exp10(x[2]))
    # j! = (J,x) -> Jac_Crit(J,eos,model,x[1]*model.params.epsilon[(1, 1)],exp10(x[2]))
    x0 = x0_crit_pure(model)
    r  = nlsolve(f!,x0)
    T_c = r.zero[1]*T̄
    v_c = exp10(r.zero[2])
    p_c = pressure(model, v_c, T_c)
    return (T_c, p_c, v_c)
end

function Obj_Crit(model::EoSModel, F, T_c, v_c)
    fun(x)  = eos(model, x, T_c,SA[1.0])
    df(x)   = ForwardDiff.derivative(fun,x)
    d2f(x)  = ForwardDiff.derivative(df,x)
    d3f(x)  = ForwardDiff.derivative(d2f,x)
    F[1] = d2f(v_c)
    F[2] = d3f(v_c)
end
## Mixture critical point solver
# function crit_mix(model::SAFT,x_c)
#     components = model.components
#     z  = create_z(model,x_c)
#     f! = (F,x) -> Obj_Crit_mix(model, F, exp10(x[2]), x[1]*prod(model.params.epsilon[i]^z[i] for i in components), x_c)
#     x0 = [1.5,log10(π/6*N_A*sum(z[i]*model.params.segment[i]*model.params.sigma[i]^3 for i in components)/0.15)]
#     r  = nlsolve(f!,x0)
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

function crit(model,z=SA[1.0])
    if length(z) == 1
        return crit_pure(model)
    else
    end
end

## Mixture saturation solver
function bubble_pressure(model, T, x; v0 =nothing)
    TYPE = promote_type(eltpype(T),eltype(x))
    
    if v0 === nothing
        y0    = 10 .*x[1,:]./(1 .+x[1,:].*(10-1))
        y0    = y0 ./sum(y0[i] for i in 1:length(x))
        X     = x[1,:]
        v0    = [log10(π/6*N_A*sum(X[i]*model.params.segment[i]*model.params.sigma[i]^3 for i in components)/0.45),
                 log10(π/6*N_A*sum(Y0[i]*model.params.segment[i]*model.params.sigma[i]^3 for i in components)/1e-4)]
        append!(v0,y0[1:end-1])
    end
    v_l   = TYPE[]
    v_v   = TYPE[]
    y     = deepcopy(x) 
    P_sat = TYPE[]
    for i in 1:first(size(x))
        f! = (F,z) -> Obj_bubble_pressure(model, F, T, exp10(z[1]), exp10(z[2]), x[i,:], z[3:end])
        j! = (J,z) -> Jac_bubble_pressure(model, J, T, exp10(z[1]), exp10(z[2]), x[i,:], z[3:end])
        r  =nlsolve(f!,j!,v0)
        append!(v_l,exp10(r.zero[1]))
        append!(v_v,exp10(r.zero[2]))
        append!(P_sat,pressure(model,v_l[i],T,x[i,:]))
        y[i,1:end-1] = r.zero[3:end]
        y[i,end] = 1-sum(r.zero[3:end])
        v0 = r.zero
    end
    return (P_sat, v_l, v_v, y)
end

function Obj_bubble_pressure(model, F, T, v_l, v_v, x, y)
    components = model.components
    append!(y,1-sum(y[i] for i in 1:(length(components)-1)))

    fun(z) = eos(model, z[end], T,z[1:end-1])
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

    fun(z) = eos(model, z[1:end-1], z[end], T)
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




#aproximates liquid volume at a known pressure and t,
#by using isothermal compressibility
function volume_compress(model,p,T,z=SA[1.0])
    v0 = vcompress_v0(model,p,T,z)
    function f_fixpoint(_v)
        _p,dpdv = p∂p∂v(model,_v,T,z)
        β = -1/_v*dpdv^-1
        _Δ =  -(p-_p)*β
        sign_Δ = sign(_Δ)
        Δ = abs(_Δ)
        vv = _v*exp(sign_Δ*Δ^(1-_Δ))#^((1+Δ)^4)
        return vv
    end
    return Solvers.fixpoint(f_fixpoint,v0,Solvers.SimpleFixPoint(),rtol = 1e-12)
    #return Roots.find_zero(f_fixpoint,v0)
end

function vcompress_v0(model,p,T,z=SA[1.0])
    lb_v   = exp10(only(lb_volume(model,z,phase=:l)))
    v0 = 1.1*lb_v
    return v0
end

function volume_virial(model,p,T, z=SA[1.] )
    B = second_virial_coeff(model,T,z)
    a = p/(R̄*T)
    b = -1
    c = -B
    return (-b + sqrt(b*b-4*a*c))/(2*a)
end