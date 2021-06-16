function try_sat_pure(model,V0,f!,T,result,error_val,converged)
    converged[] =  false
    try
        res = sat_pure(model,V0,f!,T)
        result[] = res
    catch e #normally, failures occur  near the critical point
        error_val[] = e
    end

    (P_sat,V_l,V_v) = result[]
    if abs(V_l-V_v) > 8*(eps(typeof(V_l)))
        converged[] =  true
    end

    return nothing
end


function sat_pure(model::EoSModel, T; V0 = nothing,debug=false)
    V_lb = lb_volume(model,SA[1.0])
    f! = (F,x) -> Obj_Sat(model, F, T, exp10(x[1]), exp10(x[2]),V_lb)

    if V0 === nothing
        V0 = x0_sat_pure(model,T)
    end
    nan = zero(T)/zero(T)
    res0 = (nan,nan,nan)
    result = Ref(res0)
    error_val = Ref{Any}(nothing)
    converged = Ref{Bool}(false)

    try_sat_pure(model,V0,f!,T,result,error_val,converged)
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


    T_c,P_c,V_c = crit_pure(model)
    ΔT = (T_c - T)
    ΔT <= 8*eps(ΔT) && throw(DomainError(T,"input temperature $T is too close or higher than critical temperature of the model $T_c"))
    Tr  = T/T_c
    V0 = x0_sat_pure_crit(model,T,T_c,P_c,V_c)

    try_sat_pure(model,V0,f!,T,result,error_val,converged)
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
function x0_sat_pure_crit(model,T,T_c,P_c,V_c)
    _1 = one(T)
    Tr = T/T_c
    Trm1 = _1-Tr
    Trmid = sqrt(Trm1)
    P_sat0 = (_1 -4.0*(Trm1)+4.8(Trm1*Trm1))*P_c
    #c = 1/Vr
    c_l = 1.0+2.0*Trmid + 0.4*Trm1 - 0.52*Trmid*Trm1 +0.115*Trm1*Trm1
    c_v = 1.0-2.0*Trmid + 0.4*Trm1 + 0.52*Trmid*Trm1 +0.207*Trm1*Trm1
    Vl0 = (1/c_l)*V_c
    Vv0 = (1/c_v)*V_c
    #Vl = volume_compress(model,P_sat0,T,V0=Vl0)
    #Vv = volume_compress(model,P_sat0,T,V0=Vv0)
    return [log10(Vl0),log10(Vv0)]
end

function sat_pure(model::EoSModel,V0,f!,T)
    r = Solvers.nlsolve(f!, V0,LineSearch(Newton()))
    #@show typeof(r)
    #@show f!(rand(2),r.info.zero)
    Vsol = minimizer(r)
    V_l = exp10(Vsol[1])
    V_v = exp10(Vsol[2])
    P_sat = pressure(model,V_v,T)
    return (P_sat,V_l,V_v)
end

function Obj_Sat(model::EoSModel, F, T, V_l, V_v,V_lb)
    #components = model.components
    fun(x) = eos(model, x[2], T,SA[x[1]])
    df(x)  = ForwardDiff.gradient(fun,x)
    df_l = df(SA[one(V_l*T),V_l*one(T)])
    df_v = df(SA[one(V_v),V_v*one(T)])
    (p_scale,μ_scale) = scale_sat_pure(model)
    #T̄ = T/T_scale(model)
    F[1] = (df_l[2]-df_v[2])*p_scale#*exp(5e-10*(V_l-V_v)^-2)*exp(1e-7*(V_l-V_lb)^-2)
    F[2] = (df_l[1]-df_v[1])*μ_scale#*exp(5e-10*(V_l-V_v)^-2)*exp(1e-7*(V_l-V_lb)^-2)
    return F
end
#=
function Obj_Sat(model::ABCubicModel, F, T, V_l, V_v,V_lb)
    #components = model.components
    ar(V) = a_res(model,V,T)
    function dar(V)
        a,da = Solvers.f∂f(ar,V)
        @show Vval = ForwardDiff.value(V)
        @show aval = ForwardDiff.value(da)
        @show rp = pressure(model,Vval,T)

        p = rp
        @show pval = ForwardDiff.value(p)

         Z = p*V/(R̄*T)
        ϕ = a + (Z-1) - log(Z)
        @show phival = ForwardDiff.value(ϕ)


        return p,exp(ϕ)
    end
    pl,ϕl = dar(V_l)
    pv,ϕv = dar(V_v)
    (p_scale,μ_scale) = scale_sat_pure(model)
    #T̄ = T/T_scale(model)
    F[1] = (pl-pv)*p_scale#*exp(5e-10*(V_l-V_v)^-2)*exp(1e-7*(V_l-V_lb)^-2)
    F[2] = (ϕl-ϕv)*μ_scale#*exp(5e-10*(V_l-V_v)^-2)*exp(1e-7*(V_l-V_lb)^-2)
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
    V_c = exp10(r.info.zero[2])
    p_c = pressure(model, V_c, T_c)
    return (T_c, p_c, V_c)
end

function Obj_Crit(model::EoSModel, F, T_c, V_c)
    d2p,d3p = ∂p2∂p3(model,V_c,T_c,SA[1.0])
    F[1] = -d2p
    F[2] = -d3p
    return F
end

function enthalpy_vap(model::EoSModel, T)
    (P_sat,V_l,V_v) = sat_pure(model,T)
   #= _dfl,fl =  ∂f(model,V_l,T,SA[1.0])
    _dfv,fv =  ∂f(model,V_v,T,SA[1.0])
    dVl,dTl = _dfl
    dVv,dTv = _dfv
    H_l = fl  - dVl*V_l - dTl*T
    H_v = fv  - dVv*V_v - dTv*T =#
    H_v = VT_enthalpy(model,V_v,T,z)
    H_l = VT_enthalpy(model,V_l,T,z)
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
    Vv = V_zeros(_pt,model,p,t_pred)

    V1 = first(Vv)
    #(V1 == V2) && return V1
    V2 = last(Vv)

    =#
    T0 = p_sat_t0(model,p)
    px = p

    _A(V, t) = eos(model, V, t,SA[1.0])
    V1old = zero(p)
    V2old = zero(p)
    Told = zero(p)
    V1,V2 = x0_sat_pure(model,T0)
    Tx = T0
    for i = 1:20
        if i > 1
            V1old = V1
            V2old = V2
        end
        A = z -> _A(z, Tx)
        _V1 = Threads.@spawn volume(model,$p,$Tx;phase=:liquid)
        _V2 = Threads.@spawn volume(model,$p,$Tx;phase=:gas)
        V1 = fetch(_V1)
        V2 = fetch(_V2)

        if abs(V1 - V1old) / V1 < 1e-15 && i > 1
            #println("V1 condition")
            break
        elseif abs(V2 - V2old) / V2 < 1e-15 && i > 1
            #println("V2 condition")
            break
        end

        _px(T) = (_A(V1, T) - _A(V2, T)) / (V2 - V1) - p
        #_px(T) = exp(_A(V1,T)-_A(V2,T)) - one(Tx)
        Told = Tx
        Tx = Solvers.ad_newton(_px, Tx)
        if abs(Told - Tx) < 1e-10 * Tx
            #println("P condition")
            break
        end
    end

    return (Tx,V1,V2)
end

function naive_sat_pure_p(model,p)
    T0 = p_sat_t0(model,p)
    ft0(t) = log(first(sat_pure(model,t))/p)
    tt0 =T0
    res = Solvers.ad_newton(ft0,tt0)
    T = res
    Vl = volume(model,p,T,phase=:l)
    Vv = volume(model,p,T,phase=:v)
    return (T,Vl,Vv)
end

export sat_pure, crit_pure, enthalpy_vap

function spinodals(model,T,Vx = nothing)
    T7 = (0.7)*one(T)*T_scale(model)
    #@show T/T_scale(model)
    #f0(_V) = 1/last(p∂p∂V(model,_V,T))

    lb_v = lb_volume(model)
    f0(k) = last(p∂p∂V(model,lb_v + exp(k),T))

    if Vx === nothing
    if T > T7
        sp_l7 = volume_compress(model,zero(T),T7)
        Tc,pc,Vc = crit_pure(model)
        V00 = log(sp_l7 - lb_v)
        sp_l = Roots.find_zero(f0,V00)

        sp_l = lb_v + exp(sp_l)
    else
        #zero pressure approach, the volume_compress function tends to stop at a spinodal at the liquid side
        sp_l = volume_compress(model,zero(T),T)
    end
    else
        V00 = log(Vx - lb_v)
        sp_l = Roots.find_zero(f0,V00)
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
    obtain the value at which ∂p/∂V = -0.5
    0.5 is a randomly chosen value.
    the idea is to choose the lowest value possible such
    as not pass the gas spinodal.
    with the value calculated, a pressure and then an aproximate liquid volume is calculated.

    The main equation is:
    z = 1+B/V = pV/RT
    p = RT/V + RTB/V2
    ∂p/∂V = -RT/V2 -2RTB/V3 = 0.5
    -RTV -2RTB = 0.5V3
    0 = 0.5V3 + RTV + 2RTB
    0.5V3 + 0V2 + RTV + 2RTB = 0
    (0.5,0,RT,2RTB)

    on B = B(V)
    p = RT/V + RTB/V2
    ∂p/∂V = -RT/V2 -2RT/V3*B + ∂B*RT/V2 = 0.5
    ∂p/∂V = -RT(1+∂B)/V2 -2RT/V3*B = 0.5

    -RTV -2RTB = 0.5V3
    0 = 0.5V3 + RTV + 2RTB
    0.5V3 + 0V2 + RTV + 2RTB = 0
    (0.5,0,RT,2RTB)
    =#
    #k*one(T),zero(T),R̄*T,2*B
    RT = R̄*T
    Vv = -2*B

    #Vv = Vv_vsa(B,zero(T),T,k)
    p,dpdV = p∂p∂V(model,Vv,T)
    #now, to calculate a new better aproximate of the function at that volume:
    #p = p0 + dpdV0(V-V0) / multiplying by V/RT and adding 1-1
    #pV/RT = 1-1+p0V/RT + dpdV0(V-V0)*V/RT
    #z = 1 + (-1+p0V/RT + dpdV0(V-V0)*V/RT)*V/V
    #z = 1 + B/V
    Bi(V) = (-1+p*V/RT + dpdV*(V-Vv)*V/RT)*V


    @show dpdV
    @show Bi(Vv)
    @show B
    return Vv
end
=#

#=
function pVplot(model,T)
    RT = OpenSAFT.R̄*T
    lb_v = OpenSAFT.lb_volume(model)
    @show lb_v
    b = OpenSAFT.second_virial_coefficient(model,T)
    V = range(log(1*lb_v),-2.0,length=500)
    V = exp.(V)
    peos = OpenSAFT.pressure.(model,V,T)
    pb(Vx) = RT/Vx + RT*b/(Vx*Vx)
    pvirial = pb.(V)
    V_vsa = OpenSAFT.vsa(model,T)
    spinodal_v = OpenSAFT.spinodals(model,T,V_vsa)
    lines(log10.(V),peos)
    lines!(log10.(V),pvirial)
    scatter!([log10.(V_vsa)],[pressure(model,V_vsa,T)])
    scatter!([log10.(spinodal_v)],[pressure(model,spinodal_v,T)],color=:red)
    current_figure()
end

=#
