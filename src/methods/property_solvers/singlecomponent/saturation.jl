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
    V0 = MVector((V0[1],V0[2]))
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
        V0 = MVector((V0[1],V0[2]))
        converged = try_sat_pure(model,V0,f!,T,result,error_val)   
        if converged
            return result[]
        end
    end
    #not converged, even trying with better critical aprox.
    return res0
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
    g_l,g_v = A_l - V_l*Av_l,A_v - V_v*Av_v
    (p_scale,μ_scale) = scales
    F[1] = -(Av_l-Av_v)*p_scale
    F[2] = (g_l-g_v)*μ_scale
    return F
end

function saturation_temperature(model::EoSModel,p)
    nan = zero(p)/zero(p)
    f(z) = Obj_sat_pure_T(model,z,p)
    Tc,pc,vc = crit_pure(model)
    if p>pc
        return (nan,nan,nan)
    else
        T = Roots.find_zero(f,(0.3*Tc,0.99*Tc))
        p,v_l,v_v = saturation_pressure(model,T)
        return T,v_l,v_v
    end 
end

function Obj_sat_pure_T(model,T,p)
    p̃,v_l,v_v = saturation_pressure(model,T)
    return p̃-p
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


#=
function choose_p(p1,p2)
    if p1 > 0 && p2 > 0
        p = (p1+p2)/2
    elseif p1 < 0 
        p = p2
    elseif p2 < 0
        p = p1
    elseif isnan(p1) || isnan(p2)
        p = zero(p1)/zero(p2)
    end
    return p
end
    an alternative algorithm for saturation pressure.
        
"""
    FugSucessiveSatPure(;max_iters=100,ptol= 1e-10,vtol=1e-15,V0=nothing) 

Performs a saturation pressure calculation, by sucessive substitution:
```julia
    vl = pressure(model,p,T,phase=:l)
    vv = pressure(model,p,T,phase=:v)
    p = (eos(model,vl,T) - eos(model,vv,T))/(vv-vl)
```
It's more stable than the default volume based method, but it is slower.
if no initial `V0` is passed, it will be calculated via `x0_sat_pure(model,T)`.if a number is passed, it will be used as an initial pressure.
"""
Base.@kwdef struct FugSucessiveSatPure{T}
    max_iters::Int = 100
    ptol::Float64 = 1e-10
    vtol::Float64 = 1e-15
    V0::T = nothing 
end

export FugSucessiveSatPure

function saturation_pressure(model::EoSModel,T,method::FugSucessiveSatPure)
    nan = zero(T)/zero(T)
    if isnothing(method.V0) || length(method.V0) == 2   
        if isnothing(method.V0)
            V1,V2 = x0_sat_pure(model,T)
        else
            V1,V2 = method.V0
        end
        V1,V2 = exp10(V1),exp10(V2)     
        p1,p2 = pressure(model,V1,T), pressure(model,V2,T)
    elseif (method.V0) isa Real
        p = method.V0
        p1,p2 = p
        V1 = volume(model,p,T;phase=:l)
        V2 = volume(model,p,T;phase=:v)
    else
        throw(error("invalid initial point for saturation pressure"))
    end
    T = T*one(nan)
    A(V, t) = eos(model, V, t,SA[1.0])
    p = choose_p(p1,p2)
    if isnan(p)
        T_c,P_c,V_c = crit_pure(model)
        if T_c < T
            @error "initial temperature $T greater than critical temperature $T_c. returning NaN"
        end
        V1,V2 = x0_sat_pure_crit(model,T,T_c,P_c,V_c)
        V1,V2 = exp10(V1),exp10(V2)     
        p1,p2 = pressure(model,V1,T), pressure(model,V2,T)
        p = choose_p(p1,p2)
    end
    V1old = nan
    V2old = nan
    if isnan(p)
        @error("cannot find initial pressure")
        return (nan,nan,nan)
    end
    pold = one(T)/zero(T)
    iters = method.max_iters
    ptol = method.ptol
    vtol = method.vtol
    
    function Fug(model,V,T,p)
        μres = ForwardDiff.derivative(z->eos_res(model,V,T,SA[z]),1.0)
        Z = p*V/(R̄*T)
        exp(μres/R̄/T)/Z
    end
    for i = 1:iters
        V1 = volume(model,p,T;phase=:l)
        V2 = volume(model,p,T;phase=:v)
        vvtol = max(abs(V1 - V1old)/V1, abs(V2 - V2old)/V2)
        pptol = abs(pold - p)/p
        if vvtol < vtol && i > 1
            @show i
            if check_valid_sat_pure(model,p,V1,V2,T)
                return (p,V1,V2)
            else
                return (nan,nan,nan)
            end
        end
        pold = p
        ϕl =  Fug(model,V1,T,p)
        ϕv =  Fug(model,V2,T,p)
 
        p = p*ϕl/ϕv

        if pptol < ptol && i > 1
            if check_valid_sat_pure(model,p,V1,V2,T)
                return (p,V1,V2)
            else
                return (nan,nan,nan)
            end
        end
        V1old = V1
        V2old = V2
    end
    return (nan,nan,nan)
end
=#