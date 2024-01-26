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
    vl::T
    vv::T
    crit::C
    crit_retry::Bool
    f_limit::Float64
    atol::Float64
    rtol::Float64
    max_iters::Int
end

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
    elseif !(vl === nothing) && !(vv === nothing)
        T = one(vl)/one(vv)
        vl,vv,_ = promote(vl,vv,T)
        return ChemPotVSaturation(vl,vv,crit,crit_retry,f_limit,atol,rtol,max_iters)
    else
        throw(ArgumentError("you need to specify both vl and vv."))
    end
end

ChemPotVSaturation(x::Tuple) = ChemPotVSaturation(vl = first(x),vv = last(x))
ChemPotVSaturation(x::Vector) = ChemPotVSaturation(vl = first(x),vv = last(x))

function saturation_pressure_impl(model::EoSModel, T, method::ChemPotVSaturation{Nothing,C}) where C
    TT = typeof(1.0*T*oneunit(eltype(model)))
    vl::TT,vv::TT = x0_sat_pure(model,T)
    crit = method.crit
    crit_retry = method.crit_retry
    f_limit = method.f_limit
    atol = method.atol
    rtol = method.rtol
    max_iters = method.max_iters
    new_method = ChemPotVSaturation{TT,C}(vl,vv,crit,crit_retry,f_limit,atol,rtol,max_iters)
    #@show new_method
    return saturation_pressure_impl(model,T,new_method)
end

function saturation_pressure_impl(model::EoSModel, T, method::ChemPotVSaturation{<:Number})
    vl0 = method.vl
    vv0 = method.vv
    _0 = zero(vl0*vv0*T*oneunit(eltype(model)))
    nan = _0/_0
    fail = (nan,nan,nan)
    valid_input = check_valid_2ph_input(vl0,vv0,true,T)
    if !valid_input
        return fail
    end
    ps,μs = scale_sat_pure(model)
    result,converged = try_2ph_pure_pressure(model,T,vl0,vv0,ps,μs,method)

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
    if 0.7*T_c < T < T_c
        x0 = x0_sat_pure_crit(model,T,T_c,p_c,V_c)
        vlc0,vvc0 = x0
        result,converged = try_2ph_pure_pressure(model,T,vlc0,vvc0,ps,μs,method)
        if converged
            return result
        end
    end
    #not converged, even trying with better critical aprox.
    return fail
end

export ChemPotVSaturation
