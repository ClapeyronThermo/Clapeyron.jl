function check_valid_sat_pure(model,P_sat,V_l,V_v,T)
    _,dpdvl = p∂p∂V(model,V_l,T,SA[1.0])
    _,dpdvv = p∂p∂V(model,V_v,T,SA[1.0])
    (dpdvl > 0) | (dpdvv > 0) && return false
    ε = abs(V_l-V_v)/(eps(typeof(V_l-V_v)))
    #if ΔV > ε then Vl and Vv are different values
    return ε > 5e9
end

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
    return check_valid_sat_pure(model,P_sat,V_l,V_v,T)
end

"""
    saturation_pressure(model::EoSModel, T, V0 = x0_sat_pure(model,T))

Performs a single component saturation equilibrium calculation, at the specified temperature `T`, of one mol of pure sustance specified by `model`

Returns `(p₀, Vₗ, Vᵥ)` where `p₀` is the saturation pressure (in Pa), `Vₗ` is the liquid saturation volume (in m³) and `Vᵥ` is the vapour saturation volume (in m³).

If the calculation fails, returns  `(NaN, NaN, NaN)`

`V0` is `[log10(Vₗ₀),log10(Vᵥ₀)]` , where `Vₗ₀`  and `Vᵥ₀` are initial guesses for the liquid and vapour volumes.
"""
function saturation_pressure(model::EoSModel, T, V0 = x0_sat_pure(model,T))
    !isone(length(model)) && throw(error("$model have more than one component."))
    T = T*T/T
    V01,V02 = V0
    TYPE = promote_type(typeof(T),typeof(V01),typeof(V02))
    if T isa Base.IEEEFloat # MVector does not work on non bits types, like BigFloat
        V0 = MVector((V01,V02))
    else
        V0 = SizedVector{2,typeof(first(V0))}((V01,V02))
    end
    nan = zero(TYPE)/zero(TYPE)    
    #scales = scale_sat_pure(model)
    f! = ObjSatPure(model,T)
    #f! = (F,x) -> Obj_Sat(model, F, T, exp10(x[1]), exp10(x[2]),scales)
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
        #@error "initial temperature $T greater than critical temperature $T_c. returning NaN"
    else
        x0 = x0_sat_pure_crit(model,T,T_c,p_c,V_c)
        V01,V02 = x0
        if T isa Base.IEEEFloat
            V0 = MVector((V01,V02))
        else
            V0 = SizedVector{2,typeof(first(V0))}((V01,V02))
        end
        converged = try_sat_pure(model,V0,f!,T,result,error_val)   
        if converged
            return result[]
        end
    end
    #not converged, even trying with better critical aprox.
    return res0
end

struct ObjSatPure{M,T}
    model::M
    ps::T
    mus::T
    Tsat::T
end

function ObjSatPure(model,T)
    ps,mus = scale_sat_pure(model)
    ObjSatPure(model,ps,mus,T)
end

function (f::ObjSatPure)(F,x)
    model = f.model
    scales = (f.ps,f.mus)
    T = f.Tsat
    return Obj_Sat(model, F, T, exp10(x[1]), exp10(x[2]),scales)
end
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
    return (log10(Vl0),log10(Vv0))
end
function sat_pure(model::EoSModel,V0,f!,T,method =LineSearch(Newton()))  
    r = Solvers.nlsolve(f!, V0 ,method )
    Vsol = Solvers.x_sol(r)
    V_l = exp10(Vsol[1])
    V_v = exp10(Vsol[2])
    P_sat = pressure(model,V_v,T)
    return (P_sat,V_l,V_v)
end

function Obj_Sat(model::EoSModel, F, T, V_l, V_v,scales)
    fun(_V) = eos(model, _V, T,SA[1.])
    A_l,Av_l = Solvers.f∂f(fun,V_l)
    A_v,Av_v =Solvers.f∂f(fun,V_v)
    g_l = muladd(-V_l,Av_l,A_l)
    g_v = muladd(-V_v,Av_v,A_v)
    (p_scale,μ_scale) = scales
    F[1] = -(Av_l-Av_v)*p_scale
    F[2] = (g_l-g_v)*μ_scale
    return F
end

#by default, starts right before the critical point, and descends via Clapeyron equation: (∂p/∂T)sat = ΔS/ΔV ≈ Δp/ΔT

function saturation_temperature(model::EoSModel,p,T0= nothing)
    if T0 === nothing
        T0 = x0_saturation_temperature(model,p)
    end
    TT = typeof(T0/T0)
    nan = zero(TT)/zero(TT)
    if isnan(T0)
        return (nan,nan,nan)
    end 
    cache = Ref{Tuple{TT,TT,TT,TT}}((nan,nan,nan,nan))
    f(T) = Obj_sat_pure_T(model,T,p,cache)
    T = Solvers.fixpoint(f,T0)
    _,_,v_l,v_v = cache[]
    return T,v_l,v_v
end

function x0_saturation_temperature(model,p)
    Tc,pc,vc = crit_pure(model)    
    p > 0.99999pc && (return zero(p)/zero(p))
    T99 = 0.99*Tc
    return T99
end

function Obj_sat_pure_T(model,T,p,cache)
    Told,pold,vlold,vvold = cache[]
    pii,vli,vvi = saturation_pressure(model,T)
    Δp = (p-pii)
    abs(Δp) < 4eps(p) && return T
    if Told < T
        if isnan(pii) && !isnan(pold)
            return (T+Told)/2
        end
    end
    cache[] = (T,pii,vli,vvi) 
    S_v = VT_entropy(model,vvi,T)
    S_l = VT_entropy(model,vli,T)
    ΔS = S_v - S_l
    ΔV = vvi - vli
    dpdt = ΔS/ΔV #≈ (p - pii)/(T-Tnew)
    Ti = T + Δp/dpdt
    return Ti
end

"""
    enthalpy_vap(model::EoSModel, T)

Calculates `ΔH`, the difference between saturated vapour and liquid enthalpies at temperature `T`, in J   
"""
function enthalpy_vap(model::EoSModel, T)
    (P_sat,V_l,V_v) = saturation_pressure(model,T)
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
function acentric_factor(model::EoSModel)
    T_c,p_c,_ = crit_pure(model)
    p = first(saturation_pressure(model,0.7*T_c))
    p_r = p/p_c
    return -log10(p_r) - 1.0
end



