
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



## Pure saturation conditions solver
function sat_pure(model::EoSModel, T; v0 = nothing)
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
    r = Solvers.nlsolve(f!, v0)
    #@show fieldnames(typeof(r.info))
    v_l = exp10(r.info.zero[1])
    v_v = exp10(r.info.zero[2])
    P_sat = pressure(model,exp10(r.info.zero[2]),T)
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



## Pure critical point solver
function crit_pure(model::EoSModel)
    T̄  = T_crit_pure(model)
    f! = (F,x) -> Obj_Crit(model, F, x[1]*T̄, exp10(x[2]))
    # j! = (J,x) -> Jac_Crit(J,eos,model,x[1]*model.params.epsilon[(1, 1)],exp10(x[2]))
    x0 = x0_crit_pure(model)
    r  = Solvers.nlsolve(f!,x0)
    T_c = r.info.zero[1]*T̄
    v_c = exp10(r.info.zero[2])
    p_c = pressure(model, v_c, T_c)
    return (T_c, p_c, v_c)
end

function Obj_Crit(model::EoSModel, F, T_c, v_c)
    d2p,d3p = ∂p2∂p3(model,v_c,T_c,SA[1.0])
    F[1] = -d2p
    F[2] = -d3p
    return F
end


## Mixture saturation solver
function bubble_pressure(model, T, x; v0 =nothing)
    TYPE = promote_type(eltype(T),eltype(x))
    lb_v = lb_volume(model,x)
    ts = T_scales(model,x)
    pmix = p_scale(model,x)
    #ps = p_scales(model,x)
    #Mollerup K
    #k0 =   (ps ./ p) .* exp.(5.42 .* (1.0 .- (ts ./ t)))
    k0 = 10.0
    if v0 === nothing
        y0    = k0 .*x./(1 .+x.*(k0 .- 1))
        y0    = y0 ./sum(y0)
        X     = x
        v0    = [log10(lb_v/0.45),
                 log10(lb_v/1e-4)]
        append!(v0,y0[1:end-1])
    end
    
    f! = (F,z) -> Obj_bubble_pressure(model, F, T, exp10(z[1]), exp10(z[2]), x,@view(z[3:end]),ts,pmix)
    #j! = (J,z) -> Jac_bubble_pressure(model, J, T, exp10(z[1]), exp10(z[2]), x[i,:], z[3:end])
    r  =Solvers.nlsolve(f!,v0)
    v_l = exp10(r.info.zero[1])
    v_v = exp10(r.info.zero[2])
    y = FractionVector(r.info.zero[3:end])
    P_sat = pressure(model,v_l,T,x)
    return (P_sat, v_l, v_v, y)
end

function Obj_bubble_pressure(model, F, T, v_l, v_v, x, y,ts,ps)
    #y = append!(y,1-sum(y))
    #y = [y...,1-sum(y)]
    y = FractionVector(y) #julia magic, check misc.jl
    μ_l = vt_chemical_potential(model,v_l,T,x)
    μ_v = vt_chemical_potential(model,v_v,T,y)
    p_l = pressure(model,v_l,T,x)
    p_v = pressure(model,v_v,T,y)
    
    for i in 1:length(x)
        F[i] = (μ_l[i]-μ_v[i])/(R̄*ts[i])
    end
    F[end] = (p_l-p_v)/ps
end

function enthalpy_vap(model::EoSModel, T)
    (P_sat,v_l,v_v) = sat_pure(model,T)
   #= _dfl,fl =  ∂f(model,v_l,T,SA[1.0])
    _dfv,fv =  ∂f(model,v_v,T,SA[1.0])
    dvl,dtl = _dfl
    dvv,dtv = _dfv
    H_l = fl  - dvl*v_l - dtl*T
    H_v = fv  - dvv*v_v - dtv*T =#
    H_v = vt_enthalpy(model,v_v,T,z)
    H_l = vt_enthalpy(model,v_l,T,z)
    H_vap=H_v-H_l
    return H_vap
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
    lb_v   = lb_volume(model,z,phase=:l)
    v0 = 1.1*lb_v
    return v0
end

function volume_virial(model,p,T, z=SA[1.] )
    B = second_virial_coefficient(model,T,z)
    a = p/(R̄*T)
    b = -1
    c = -B
    Δ = b*b-4*a*c
    n = sum(z)
    if Δ <= 0
        #virial approximation could not be calculated
        #degrade to ideal aprox
        return n*R̄*T/p 
        
    end 
    return (-b + sqrt(b*b-4*a*c))/(2*a)
end

#(z = pv/rt)
#(RT/p = v/z)