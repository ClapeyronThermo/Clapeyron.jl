function sat_pure(model::EoSModel, T; v0 = nothing)
    v_lb = lb_volume(model,SA[1.0])
    f! = (F,x) -> Obj_Sat(model, F, T, exp10(x[1]), exp10(x[2]),v_lb)
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
        (P_sat,v_l,v_v) = sat_pure(model,v0,f!,T)
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

function Obj_Sat(model::EoSModel, F, T, v_l, v_v,v_lb)
    #components = model.components
    fun(x) = eos(model, x[2], T,SA[x[1]])
    df(x)  = ForwardDiff.gradient(fun,x)
    df_l = df(SA[one(v_l*T),v_l])
    df_v = df(SA[one(v_v*T),v_v])
    (p_scale,μ_scale) = scale_sat_pure(model)
    #T̄ = T/T_scale(model)
    F[1] = (df_l[2]-df_v[2])*p_scale*exp(5e-10*(v_l-v_v)^-2)*exp(1e-7*(v_l-v_lb)^-2)
    F[2] = (df_l[1]-df_v[1])*μ_scale*exp(5e-10*(v_l-v_v)^-2)*exp(1e-7*(v_l-v_lb)^-2)
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
<<<<<<< HEAD
        sp_l7 = volume_compress(model,zero(T),T7)
        f0(k) = last(p∂p∂v(model,lb_v + exp(k),T))

=======
        sp_l7 = volume_compress(model,zero(T),T7)
>>>>>>> 43be74673998cbf47e35d31ce01100ff65d8d76e
        tc,pc,vc = crit_pure(model)
        v00 = log(sp_l7 - lb_v)
        sp_l = Roots.find_zero(f0,v00)

        sp_l = lb_v + exp(sp_l)
    else
        #zero pressure approach, the volume_compress function tends to stop at a spinodal at the liquid side
        sp_l = volume_compress(model,zero(T),T)
    end
<<<<<<< HEAD
    return (sp_l)#sp_v)
end
=======
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
function vsa(model,T)
    k = -0.5*one(T)#*(p_scale(model)/lb_volume(model))
    B = second_virial_coefficient(model,T)
    #=basis
    obtain the value at which ∂p/∂v = -0.5
    0.5 is a randomly chosen value.
    the idea is to choose the lowest value possible such
    as not pass the gas spinodal.
    with the value calculated, a pressure and then an aproximate liquid volume is calculated.

    The main equation is:
    z = 1+B/v = pv/RT
    p = RT/v + RTB/v2
    ∂p/∂v = -RT/v2 -2RTB/v3 = 0.5
    -RTv -2RTB = 0.5v3
    0 = 0.5v3 + RTv + 2RTB
    0.5v3 + 0v2 + RTv + 2RTB = 0
    (0.5,0,RT,2RTB)

    on B = B(v)
    p = RT/v + RTB/v2
    ∂p/∂v = -RT/v2 -2RT/v3*B + ∂B*RT/v2 = 0.5
    ∂p/∂v = -RT(1+∂B)/v2 -2RT/v3*B = 0.5

    -RTv -2RTB = 0.5v3
    0 = 0.5v3 + RTv + 2RTB
    0.5v3 + 0v2 + RTv + 2RTB = 0
    (0.5,0,RT,2RTB)
    =#
    #k*one(T),zero(T),R̄*T,2*B
    RT = R̄*T
    vv = -2*B

    #vv = vv_vsa(B,zero(T),T,k)
    p,dpdv = p∂p∂v(model,vv,T)
    #now, to calculate a new better aproximate of the function at that volume:
    #p = p0 + dpdv0(v-v0) / multiplying by v/RT and adding 1-1
    #pv/RT = 1-1+p0v/RT + dpdv0(v-v0)*v/RT
    #z = 1 + (-1+p0v/RT + dpdv0(v-v0)*v/RT)*v/v
    #z = 1 + B/v
    Bi(v) = (-1+p*v/RT + dpdv*(v-vv)*v/RT)*v


    @show dpdv
    @show Bi(vv)
    @show B
    return vv
end

function vv_vsa(B,db,T,k)
    RT = R̄*T*(one(T)+db)
    resv = Solvers.roots3(2*RT*B,RT,zero(T),one(T)*k)
    #if k =0
    #RTv + 2RTB = 0
    #v + 2B = 0
    #v = -2B

    #@show real.(resv)
    vv = sort(real.(resv))[2]
    return vv
end
function vsa7(model)
    T = 0.7*T_scale(model)
    return vsa(model,T)
end

function vsak(model,k)
    T = k*T_scale(model)
    return vsa(model,T)
end

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
>>>>>>> 43be74673998cbf47e35d31ce01100ff65d8d76e
