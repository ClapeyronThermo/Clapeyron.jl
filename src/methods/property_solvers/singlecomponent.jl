function try_sat_pure(model,V0,f!,T,result,error_val,method = LineSearch(Newton()))
    if !isfinite(V0[1]) | !isfinite(V0[2])
        return false
    end
    try
        res = sat_pure(model,V0,f!,T,method)
        result[] = res
    catch e #normally, failures occur near the critical point
        error_val[] = e
        return false
    end

    (P_sat,V_l,V_v) = result[]
    _,dpdvl = p∂p∂V(model,V_l,T,SA[1.0])
    _,dpdvv = p∂p∂V(model,V_v,T,SA[1.0])
    (dpdvl > 0) | (dpdvv > 0) && return false
    ε = abs(V_l-V_v)/(eps(typeof(V_l-V_v)))
    ups = (eps(typeof(V_l-V_v)))
    #@show abs(V_l-V_v),ups
    if (ε > 5e9)
        #if ΔV > ε then Vl and Vv are different values
        return true
    end

    return false
end

"""
    sat_pure(model::EoSModel, T, V0 = x0_sat_pure(model,T))

Performs a single component saturation equilibrium calculation, at the specified temperature `T`, of one mol of pure sustance specified by `model`

Returns `(p₀, Vₗ, Vᵥ)` where `p₀` is the saturation pressure (in Pa), `Vₗ` is the liquid saturation volume (in m³) and `Vᵥ` is the vapour saturation volume (in m³).

If the calculation fails, returns  `(NaN, NaN, NaN)`

`V0` is `[log10(Vₗ₀),log10(Vᵥ₀)]` , where `Vₗ₀`  and `Vᵥ₀` are initial guesses for the liquid and vapour volumes.
"""
function sat_pure(model::EoSModel, T, V0 = x0_sat_pure(model,T))
    T = T*T/T
    V_lb = lb_volume(model,SA[1.0])
    TYPE = promote_type(typeof(T),typeof(V_lb))
    nan = zero(TYPE)/zero(TYPE)    
    scales = scale_sat_pure(model)
    f! = (F,x) -> Obj_Sat(model, F, T, exp10(x[1]), exp10(x[2]),scales)
    res0 = (nan,nan,nan)
    result = Ref(res0)
    error_val = Ref{Any}(nothing)
    converged = try_sat_pure(model,V0,f!,T,result,error_val)   
    #did not converge, but didnt error.
    if converged
        return result[]
    end
    (T_c, p_c, V_c) = crit_pure(model)
    if abs(T_c-T) < eps(typeof(T))
        return (p_c,V_c,V_c)
    end
    if T_c < T
        @error "initial temperature $T greater than critical temperature $T_c. returning NaN"
    else
        V0 = x0_sat_pure_crit(model,T,T_c,p_c,V_c)
        converged = try_sat_pure(model,V0,f!,T,result,error_val)   
        if converged
            return result[]
        end
    end
    #not converged, even trying with better critical aprox.
    return res0
end


 #   x0_sat_pure_crit(model,T,T_c,P_c,V_c)



#with the critical point, we can perform a
#corresponding states approximation with the
#propane reference equation of state
function x0_sat_pure_crit(model,T,T_c,P_c,V_c) 
    h = V_c*5000
    T0 = 369.89*T/T_c
    Vl0 = (1.0/_propaneref_rholsat(T0))*h
    Vv0 = (1.0/_propaneref_rhovsat(T0))*h
    _1 = SA[1.0]
    #μ_l = only(VT_chemical_potential(model,Vl0,T,_1))
    #μ_v = only(VT_chemical_potential(model,Vv0,T,_1))
    #@show (μ_l < μ_v,T/T_c)
    #if μ_l < μ_v 
      #@show μ_l,μ_v
    #end
    # _,dpdvv = p∂p∂V(model,Vv0,T,SA[1.0])
    # @show dpdvv*Vv0
    # _,dpdvv = p∂p∂V(model,2*Vv0,T,SA[1.0])
    # @show dpdvv*Vv0
    return [log10(Vl0),log10(Vv0)]
end
function sat_pure(model::EoSModel,V0,f!,T,method =LineSearch(Newton()) )  
    r = Solvers.nlsolve(f!, V0 ,method )
    Vsol = Solvers.x_sol(r)
    V_l = exp10(Vsol[1])
    V_v = exp10(Vsol[2])
    P_sat = pressure(model,V_v,T)
    return (P_sat,V_l,V_v)
end

function Obj_Sat(model::EoSModel, F, T, V_l, V_v,scales)
    fun(x) = eos(model, x[2], T,SA[x[1]])
    df(x)  = ForwardDiff.gradient(fun,x)
    df_l = df(SA[one(V_l*T),V_l*one(T)])
    df_v = df(SA[one(V_v),V_v*one(T)])
    (p_scale,μ_scale) = scales
    #T̄ = T/T_scale(model)
    F[1] = (df_l[2]-df_v[2])*p_scale
    F[2] = (df_l[1]-df_v[1])*μ_scale
    return F
   
    #=
    This is the original algorithm, but it seems
    that evaluating both pressure and chemical potential 
    at the same time is faster.
    
    #GERG2008(["water"]) - ForwardDiff one pass - 74.100 μs
    #GERG2008(["water"]) - original - 108.000 μs
    _1 = SA[1.0]
    μ_l = only(VT_chemical_potential(model,V_l,T,_1))
    μ_v = only(VT_chemical_potential(model,V_v,T,_1))
    p_l = -pressure(model,V_l,T,_1)
    p_v = -pressure(model,V_v,T,_1)
    F[1] = (μ_l-μ_v)*p_scale
    F[2] = (p_l-p_v)*μ_scale
    return F
    =#
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


"""
    crit_pure(model::EoSModel,x0=nothing)

Calculates the critical point of a single component modelled by `model`. 

Returns `(Tc, pc, Vc)` where `Tc` is the critical temperature (in K), `pc` is the critical pressure (in Pa) and `Vc` is the critical volume (in  m³)
"""
function crit_pure(model::EoSModel,x0=nothing)
    T̄  = T_scale(model)
    f! = (F,x) -> obj_crit(model, F, x[1]*T̄, exp10(x[2]))
    if x0 === nothing
        x0 = x0_crit_pure(model)
    end

    solver_res = Solvers.nlsolve(f!, x0)
    #print(solver_res)
    r  = Solvers.x_sol(solver_res)
    T_c = r[1]*T̄
    V_c = exp10(r[2])
    p_c = pressure(model, V_c, T_c)
    return (T_c, p_c, V_c)
end

function obj_crit(model::EoSModel, F, T_c, V_c)
    _1 = one(T_c+V_c)
    ∂²A∂V², ∂³A∂V³ = ∂²³f(model, V_c, T_c, SA[_1])
    F[1] = -∂²A∂V²
    F[2] = -∂³A∂V³
    return F
end
"""
    enthalpy_vap(model::EoSModel, T)

Calculates `ΔH`, the difference between saturated vapour and liquid enthalpies at temperature `T`, in J   
"""
function enthalpy_vap(model::EoSModel, T)
    (P_sat,V_l,V_v) = sat_pure(model,T)
    H_v = VT_enthalpy(model,V_v,T)
    H_l = VT_enthalpy(model,V_l,T)
    H_vap=H_v -H_l
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


# function p_sat_t0(model,p)
#     #ln(pr) = h*(1-1/tr)
#     #1/(1-ln(pr)/h)*tc = t
#     tc = T_scale(model)
#     pc = p_scale(model)
#     scale = 0.7
#     t7 = scale*tc
#     p7 = first(sat_pure(model,t7))
#     h = 2.3333333333333335*log(pc/p7)
#     #T/Tc = 1(1-log(p/pc)/h)
#     return 1/(1-log(p/pc)/h)*tc
# end

# function sat_pure_p(model,p)
#     #=
#     model_pred = single_sat_aprox(model)
#     t_pred = p_sat_t0(model,p) #prediction temperature
#     Vv = V_zeros(_pt,model,p,t_pred)

#     V1 = first(Vv)
#     #(V1 == V2) && return V1
#     V2 = last(Vv)

#     =#
#     T0 = p_sat_t0(model,p)
#     px = p

#     _A(V, t) = eos(model, V, t,SA[1.0])
#     V1old = zero(p)
#     V2old = zero(p)
#     Told = zero(p)
#     V1,V2 = x0_sat_pure(model,T0)
#     Tx = T0
#     for i = 1:20
#         if i > 1
#             V1old = V1
#             V2old = V2
#         end
#         A = z -> _A(z, Tx)
#         _V1 = Threads.@spawn volume(model,$p,$Tx;phase=:liquid)
#         _V2 = Threads.@spawn volume(model,$p,$Tx;phase=:gas)
#         V1 = fetch(_V1)
#         V2 = fetch(_V2)

#         if abs(V1 - V1old) / V1 < 1e-15 && i > 1
#             #println("V1 condition")
#             break
#         elseif abs(V2 - V2old) / V2 < 1e-15 && i > 1
#             #println("V2 condition")
#             break
#         end

#         _px(T) = (_A(V1, T) - _A(V2, T)) / (V2 - V1) - p
#         #_px(T) = exp(_A(V1,T)-_A(V2,T)) - one(Tx)
#         Told = Tx
#         Tx = Solvers.ad_newton(_px, Tx)
#         if abs(Told - Tx) < 1e-10 * Tx
#             #println("P condition")
#             break
#         end
#     end

#     return (Tx,V1,V2)
# end

# function naive_sat_pure_p(model,p)
#     T0 = p_sat_t0(model,p)
#     ft0(t) = log(first(sat_pure(model,t))/p)
#     tt0 =T0
#     res = Solvers.ad_newton(ft0,tt0)
#     T = res
#     Vl = volume(model,p,T,phase=:l)
#     Vv = volume(model,p,T,phase=:v)
#     return (T,Vl,Vv)
# end

export sat_pure, crit_pure, enthalpy_vap

# function spinodals(model,T,Vx = nothing)
#     T7 = (0.7)*one(T)*T_scale(model)
#     #@show T/T_scale(model)
#     #f0(_V) = 1/last(p∂p∂V(model,_V,T))

#     lb_v = lb_volume(model)
#     f0(k) = last(p∂p∂V(model,lb_v + exp(k),T))

#     if Vx === nothing
#     if T > T7
#         sp_l7 = volume_compress(model,zero(T),T7)
#         Tc,pc,Vc = crit_pure(model)
#         V00 = log(sp_l7 - lb_v)
#         sp_l = Roots.find_zero(f0,V00)

#         sp_l = lb_v + exp(sp_l)
#     else
#         #zero pressure approach, the volume_compress function tends to stop at a spinodal at the liquid side
#         sp_l = volume_compress(model,zero(T),T)
#     end
#     else
#         V00 = log(Vx - lb_v)
#         sp_l = Roots.find_zero(f0,V00)
#         sp_l = lb_v + exp(sp_l)

#     end
#     #now calculate gas spinodal
#     #@show Pmin = pressure(model,sp_l,T)
#     #p could be negative, and a volume search fails here.
#     #@show volume_virial(model,Pmin,T)
#     return sp_l#sp_v)
# end

