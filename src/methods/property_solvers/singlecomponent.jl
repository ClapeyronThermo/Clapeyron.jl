function try_sat_pure(model,v0,f!,T,result,error_val,converged)
    converged[] =  false
    try
        res = sat_pure(model,v0,f!,T)
        result[] = res
    catch e #normally, failures occur  near the critical point
        error_val[] = e
    end

    (P_sat,v_l,v_v) = result[]
    if abs(v_l-v_v) > 8*(eps(typeof(v_l)))
        converged[] =  true
    end        

    return nothing
end


function sat_pure(model::EoSModel, T; v0 = nothing,debug=false)
    v_lb = lb_volume(model,SA[1.0])
    f! = (F,x) -> Obj_Sat(model, F, T, exp10(x[1]), exp10(x[2]),v_lb)

    if v0 === nothing
        v0 = x0_sat_pure(model,T)
    end
    nan = zero(T)/zero(T)
    res0 = (nan,nan,nan)
    result = Ref(res0)
    error_val = Ref{Any}(nothing)
    converged = Ref{Bool}(false)

    try_sat_pure(model,v0,f!,T,result,error_val,converged)
    if converged[]
        return result[]
    end
    if debug
        if error_val[] !== nothing
            @warn "initial saturation calculation failed with error $error_val[]"
        end
    end
    result[] = res0
    #convergence not achieved, trying critical aproximation
   
    
    T_c,P_c,v_c = crit_pure(model)
    ΔT = (T_c - T)
    ΔT <= 8*eps(ΔT) && throw(DomainError(T,"input temperature $T is too close or higher than critical temperature of the model $T_c"))
    Tr  = T/T_c
    v0 = x0_sat_pure_crit(model,T,T_c,P_c,v_c)
    
    try_sat_pure(model,v0,f!,T,result,error_val,converged)
    if converged[]
        return result[]
    else
        @warn "the procedure converged to a trivial value at T=$T"
        return result[]
    end

    
    if debug
        throw(error_val[])
    else
        throw("unable to calculate equilibria at T=$T")
    end
    

     
end


#=
based on:
DOI: 10.1007/s10910-007-9272-4
Journal of Mathematical Chemistry, Vol. 43, No. 4, May 2008 (© 2007)
The van der Waals equation: analytical and approximate solutions
=#
function x0_sat_pure_crit(model,T,T_c,P_c,v_c)
    _1 = one(T)
    Tr = T/T_c
    Trm1 = _1-Tr
    Trmid = sqrt(Trm1)
    P_sat0 = (_1 -4.0*(Trm1)+4.8(Trm1*Trm1))*P_c
    #c = 1/vr
    c_l = 1.0+2.0*Trmid + 0.4*Trm1 - 0.52*Trmid*Trm1 +0.115*Trm1*Trm1
    c_v = 1.0-2.0*Trmid + 0.4*Trm1 + 0.52*Trmid*Trm1 +0.207*Trm1*Trm1
    vl0 = (1/c_l)*v_c
    vv0 = (1/c_v)*v_c
    #vl = volume_compress(model,P_sat0,T,v0=vl0)
    #vv = volume_compress(model,P_sat0,T,v0=vv0)
    return [log10(vl0),log10(vv0)]
end

function sat_pure(model::EoSModel,v0,f!,T)
    r = Solvers.nlsolve(f!, v0,LineSearch(Newton()))
    #@show typeof(r)
    #@show f!(rand(2),r.info.zero)
    vsol = minimizer(r)
    v_l = exp10(vsol[1])
    v_v = exp10(vsol[2])
    P_sat = pressure(model,v_v,T)
    return (P_sat,v_l,v_v)
end

function Obj_Sat(model::EoSModel, F, T, v_l, v_v,v_lb)
    #components = model.components
    fun(x) = eos(model, x[2], T,SA[x[1]])
    df(x)  = ForwardDiff.gradient(fun,x)
    df_l = df(SA[one(v_l*T),v_l*one(T)])
    df_v = df(SA[one(v_v),v_v*one(T)])
    (p_scale,μ_scale) = scale_sat_pure(model)
    #T̄ = T/T_scale(model)
    F[1] = (df_l[2]-df_v[2])*p_scale#*exp(5e-10*(v_l-v_v)^-2)*exp(1e-7*(v_l-v_lb)^-2)
    F[2] = (df_l[1]-df_v[1])*μ_scale#*exp(5e-10*(v_l-v_v)^-2)*exp(1e-7*(v_l-v_lb)^-2)
    return F
end
#=
function Obj_Sat(model::ABCubicModel, F, T, v_l, v_v,v_lb)
    #components = model.components
    ar(v) = a_res(model,v,T)
    function dar(v) 
        a,da = Solvers.f∂f(ar,v)
        @show vval = ForwardDiff.value(v)
        @show aval = ForwardDiff.value(da)
        @show rp = pressure(model,vval,T)

        p = rp
        @show pval = ForwardDiff.value(p)

         Z = p*v/(R̄*T)
        ϕ = a + (Z-1) - log(Z)
        @show phival = ForwardDiff.value(ϕ)


        return p,exp(ϕ)
    end
    pl,ϕl = dar(v_l)
    pv,ϕv = dar(v_v)
    (p_scale,μ_scale) = scale_sat_pure(model)
    #T̄ = T/T_scale(model)
    F[1] = (pl-pv)*p_scale#*exp(5e-10*(v_l-v_v)^-2)*exp(1e-7*(v_l-v_lb)^-2)
    F[2] = (ϕl-ϕv)*μ_scale#*exp(5e-10*(v_l-v_v)^-2)*exp(1e-7*(v_l-v_lb)^-2)
    return F
end
=#


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



"""
    acentric_factor(model::EoSModel)

calculates the acentric factor using its definition:

    ω = -log10(psatᵣ) -1, at Tᵣ = 0.7
To do so, it calculates the critical temperature (using `crit_pure`) and performs a saturation calculation (with `sat_pure`)

"""
function acentric_factor(model)
    T_c,p_c,_ = crit_pure(model)
    T = 0.7*T_c
    p = first(sat_pure(model,T))
    p_r = p/p_c
    return -log10(p_r) - 1.0
end


function p_sat_t0(model,p)
    #ln(pr) = h*(1-1/tr)
    #1/(1-ln(pr)/h)*tc = t
    tc = T_scale(model)
    pc = p_scale(model)
    scale = 0.7
    t7 = scale*tc
    p7 = first(sat_pure(model,t7))
    h = 2.3333333333333335*log(pc/p7)
    #T/Tc = 1(1-log(p/pc)/h)
    return 1/(1-log(p/pc)/h)*tc
end

function sat_pure_p(model,p)
    #=
    model_pred = single_sat_aprox(model)
    t_pred = p_sat_t0(model,p) #prediction temperature
    vv = v_zeros(_pt,model,p,t_pred)

    v1 = first(vv)
    #(v1 == v2) && return v1
    v2 = last(vv)

    =#
    T0 = p_sat_t0(model,p)
    px = p

    _A(v, t) = eos(model, v, t,SA[1.0])
    v1old = zero(p)
    v2old = zero(p)
    Told = zero(p)
    v1,v2 = x0_sat_pure(model,T0)
    Tx = T0
    for i = 1:20
        if i > 1
            v1old = v1
            v2old = v2
        end
        A = z -> _A(z, Tx)
        _v1 = Threads.@spawn volume(model,$p,$Tx;phase=:liquid)
        _v2 = Threads.@spawn volume(model,$p,$Tx;phase=:gas)
        v1 = fetch(_v1)
        v2 = fetch(_v2)

        if abs(v1 - v1old) / v1 < 1e-15 && i > 1
            #println("v1 condition")
            break
        elseif abs(v2 - v2old) / v2 < 1e-15 && i > 1
            #println("v2 condition")
            break
        end

        _px(T) = (_A(v1, T) - _A(v2, T)) / (v2 - v1) - p
        #_px(T) = exp(_A(v1,T)-_A(v2,T)) - one(Tx)
        Told = Tx
        Tx = Solvers.ad_newton(_px, Tx)
        if abs(Told - Tx) < 1e-10 * Tx
            #println("P condition")
            break
        end
    end

    return (Tx,v1,v2)
end

function naive_sat_pure_p(model,p)
    T0 = p_sat_t0(model,p)
    ft0(t) = log(first(sat_pure(model,t))/p)
    tt0 =T0
    res = Solvers.ad_newton(ft0,tt0)
    T = res
    vl = volume(model,p,T,phase=:l)
    vv = volume(model,p,T,phase=:v)
    return (T,vl,vv)
end

export sat_pure, crit_pure, enthalpy_vap

function spinodals(model,T,vx = nothing)
    T7 = (0.7)*one(T)*T_scale(model)
    #@show T/T_scale(model)
    #f0(_v) = 1/last(p∂p∂v(model,_v,T))
    
    lb_v = lb_volume(model)
    f0(k) = last(p∂p∂v(model,lb_v + exp(k),T))

    if vx === nothing
    if T > T7
        sp_l7 = volume_compress(model,zero(T),T7) 
        tc,pc,vc = crit_pure(model)
        v00 = log(sp_l7 - lb_v)
        sp_l = Roots.find_zero(f0,v00) 

        sp_l = lb_v + exp(sp_l)
    else
        #zero pressure approach, the volume_compress function tends to stop at a spinodal at the liquid side
        sp_l = volume_compress(model,zero(T),T)
    end
    else
        v00 = log(vx - lb_v)
        sp_l = Roots.find_zero(f0,v00) 
        sp_l = lb_v + exp(sp_l)

    end
    #now calculate gas spinodal
    #@show Pmin = pressure(model,sp_l,T)
    #p could be negative, and a volume search fails here.
    #@show volume_virial(model,Pmin,T)
    return sp_l#sp_v)
end


#=
function pvplot(model,T)
    RT = OpenSAFT.R̄*T
    lb_v = OpenSAFT.lb_volume(model)
    @show lb_v
    b = OpenSAFT.second_virial_coefficient(model,T)
    v = range(log(1*lb_v),-2.0,length=500)
    v = exp.(v)
    peos = OpenSAFT.pressure.(model,v,T)
    pb(vx) = RT/vx + RT*b/(vx*vx)
    pvirial = pb.(v)
    v_vsa = OpenSAFT.vsa(model,T)
    spinodal_v = OpenSAFT.spinodals(model,T,v_vsa)
    lines(log10.(v),peos)
    lines!(log10.(v),pvirial)
    scatter!([log10.(v_vsa)],[pressure(model,v_vsa,T)])
    scatter!([log10.(spinodal_v)],[pressure(model,spinodal_v,T)],color=:red)
    current_figure()
end

=#
