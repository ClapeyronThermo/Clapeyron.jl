#TODO: better name

struct ChemPotVSaturation{T} <: SaturationMethod
    vl::Union{Nothing,T}
    vv::Union{Nothing,T}
end

ChemPotVSaturation(x::Tuple) = ChemPotVSaturation(first(x),last(x))
ChemPotVSaturation(x::Vector) = ChemPotVSaturation(first(x),last(x))

function vec2(method::ChemPotVSaturation{T},opt = true) where T <:Real
    return vec2(method.vl,method.vv,t)
end

function ChemPotVSaturation(;vl = nothing,vv = nothing)
    if (vl === nothing) && (vv === nothing)
        return ChemPotVSaturation{Nothing}(nothing,nothing)
    elseif !(vl === nothing) && (vv === nothing)
        vl = float(vl)
        return ChemPotVSaturation(vl,vv)
    elseif (vl === nothing) && !(vv === nothing)
        vv = float(vv)
        return ChemPotVSaturation(vl,vv)
    else
        T = one(vl)/one(vv)
        vl,vv,_ = promote(vl,vv,T)
        return ChemPotVSaturation(vl,vv)
    end
end

function check_valid_sat_pure(model,P_sat,V_l,V_v,T)
    _,dpdvl = p∂p∂V(model,V_l,T,SA[1.0])
    _,dpdvv = p∂p∂V(model,V_v,T,SA[1.0])
    (dpdvl > 0) | (dpdvv > 0) && return false
    ε = abs(V_l-V_v)/(eps(typeof(V_l-V_v)))
    #if ΔV > ε then Vl and Vv are different values
    return ε > 5e7
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

function saturation_pressure(model,T,V0::Union{Tuple,Vector} = x0_sat_pure(model,T))
    method = ChemPotVSaturation(V0)
    return saturation_pressure_impl(model,T,method)
end

function saturation_pressure_impl(model::EoSModel, T, method::ChemPotVSaturation)
    V0 = vec2(method,T)
    V01,V02 = V0
    TYPE = eltype(V0)
    nan = zero(TYPE)/zero(TYPE)    
    f! = ObjSatPure(model,T) #functor
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
