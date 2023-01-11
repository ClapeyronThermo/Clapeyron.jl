#TODO: better name
"""
    ChemPotVSaturation <: SaturationMethod
    ChemPotVSaturation(V0)
    ChemPotVSaturation(;vl = nothing,
                        vv = nothing,
                        crit = nothing,
                        crit_retry = true
                        f_limit = 0.0,
                        atol = 1e-8,
                        rtol = 1e-12,
                        max_iters = 10^4)

Default `saturation_pressure` Saturation method used by `Clapeyron.jl`. It uses equality of Chemical Potentials with a volume basis. If no volumes are provided, it will use  [`x0_sat_pure`](@ref).

If those initial guesses fail and the specification is near critical point, it will try one more time, using Corresponding States instead.

when `crit_retry` is true, if the initial solve fail, it will try to obtain a better estimate by calculating the critical point. 

`f_limit`, `atol`, `rtol`, `max_iters` are passed to the non linear system solver.
"""
struct ChemPotVSaturation{T,C} <: SaturationMethod
    vl::Union{Nothing,T}
    vv::Union{Nothing,T}
    crit::C
    crit_retry::Bool
    f_limit::Float64
    atol::Float64
    rtol::Float64
    max_iters::Int
end

ChemPotVSaturation(x::Tuple) = ChemPotVSaturation(vl = first(x),vv = last(x))
ChemPotVSaturation(x::Vector) = ChemPotVSaturation(vl = first(x),vv = last(x))

function ChemPotVSaturation(;vl = nothing,
                            vv = nothing,
                            crit = nothing,
                            crit_retry = true,
                            f_limit = 0.0,
                            atol = 1e-8,
                            rtol = 1e-12,
                            max_iters = 10000)

    if (vl === nothing) && (vv === nothing)
        return ChemPotVSaturation{Nothing,typeof(crit)}(nothing,nothing,crit,crit_retry,f_limit,atol,rtol,max_iters)
    elseif !(vl === nothing) && (vv === nothing)
        vl = float(vl)
        return ChemPotVSaturation(vl,vv,crit,crit_retry,f_limit,atol,rtol,max_iters)
    elseif (vl === nothing) && !(vv === nothing)
        vv = float(vv)
        return ChemPotVSaturation(vl,vv,crit,crit_retry,f_limit,atol,rtol,max_iters)
    else
        T = one(vl)/one(vv)
        vl,vv,_ = promote(vl,vv,T)
        return ChemPotVSaturation(vl,vv,crit,crit_retry,f_limit,atol,rtol,max_iters)
    end
end

function saturation_pressure(model::EoSModel,T,V0::Union{Tuple,Vector})
    single_component_check(saturation_pressure,model)
    method = ChemPotVSaturation(V0)
    T = T*T/T
    return saturation_pressure_impl(model,T,method)
end

function saturation_pressure(model::EoSModel,T)
   saturation_pressure(model,T,ChemPotVSaturation())
end

function saturation_pressure_impl(model::EoSModel, T, method::ChemPotVSaturation{Nothing})
    vl,vv = x0_sat_pure(model,T)
    crit = method.crit
    return saturation_pressure_impl(model,T,ChemPotVSaturation(;vl,vv,crit))
end

function saturation_pressure_impl(model::EoSModel, T, method::ChemPotVSaturation{<:Number})
    V0 = vec2(log(method.vl),log(method.vv),T)
    V01,V02 = V0
    TYPE = eltype(V0)
    nan = zero(TYPE)/zero(TYPE)
    f! = ObjSatPure(model,T) #functor
    fail = (nan,nan,nan)
    if !isfinite(V0[1]) | !isfinite(V0[2]) | !isfinite(T)
        #error in initial conditions
        return fail
    end
    result,converged = sat_pure(f!,V0,method)
    #did not converge, but didnt error.
    if converged
        return result
    end
    if !method.crit_retry
        return fail
    end

    crit = method.crit
    if isnothing(crit)
        crit = crit_pure(model)
    end
    T_c, p_c, V_c = crit
    if abs(T_c-T) < eps(typeof(T))
        return (p_c,V_c,V_c)
    end
        #@error "initial temperature $T greater than critical temperature $T_c. returning NaN"
    if T < T_c
        x0 = x0_sat_pure_crit(model,T,T_c,p_c,V_c)
        V01,V02 = x0
        V0 = vec2(log(V01),log(V02),T)
        result,converged = sat_pure(f!,V0,method)
        if converged
            return result
        end
    end
    #not converged, even trying with better critical aprox.
    return fail
end

struct ObjSatPure{M,T}
    model::M
    ps::T
    mus::T
    Tsat::T
end

function ObjSatPure(model,T)
    ps,mus = scale_sat_pure(model)
    ps,mus,T = promote(ps,mus,T)
    ObjSatPure(model,ps,mus,T)
end

function (f::ObjSatPure)(F,x)
    model = f.model
    scales = (f.ps,f.mus)
    T = f.Tsat
    return Obj_Sat(model, F, T, exp(x[1]), exp(x[2]),scales)
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
    return Vl0,Vv0
end

function sat_pure(model,T,V0,method)
    f! = ObjSatPure(model,T)
    return sat_pure(f!,V0,method)
end

function sat_pure(f!::ObjSatPure,V0,method)
    model, T = f!.model, f!.Tsat
    nan = zero(eltype(V0))/zero(eltype(V0))
    if !isfinite(V0[1]) | !isfinite(V0[2]) | !isfinite(T)
        return (nan,nan,nan), false
    end
    r = Solvers.nlsolve(f!, V0 ,LineSearch(Newton()),NEqOptions(method),ForwardDiff.Chunk{2}())
    Vsol = Solvers.x_sol(r)
    V_l = exp(Vsol[1])
    V_v = exp(Vsol[2])
    P_sat = pressure(model,V_v,T)
    res = (P_sat,V_l,V_v)
    valid = check_valid_sat_pure(model,P_sat,V_l,V_v,T)
    return res,valid
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

export ChemPotVSaturation
