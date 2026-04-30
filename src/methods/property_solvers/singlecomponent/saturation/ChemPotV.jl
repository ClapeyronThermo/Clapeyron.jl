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
Default `saturation_pressure` Saturation method used by `Clapeyron.jl`.
It uses equality of Chemical Potentials with a volume basis.
If no volumes are provided, it will use  [`x0_sat_pure`](@ref).
When `crit_retry` is true, if the initial solve fail, it will try to obtain a better estimate by calculating the critical point.
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
    TT = Base.promote_eltype(model,T)
    crit = method.crit
    if crit !== nothing && !has_fast_crit_pure(model)
        vl,vv = x0_sat_pure(model,T,crit)
    else
        vl,vv = x0_sat_pure(model,T)
    end
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
    valid_input = _is_positive((vl0,vv0,T))
    if !valid_input
        return fail
    end
    ps,μs = equilibria_scale(model)
    result,converged = try_2ph_pure_pressure(model,T,vl0,vv0,ps,μs,method)

    crit = method.crit
    converged && return result #converged result.
    !converged && crit !== nothing && return fail #we already used critical information
    !method.crit_retry && return fail #we dont want to try again
    if isnothing(crit)
        crit = crit_pure(model)
    end
    T_c, p_c, V_c = crit
    if abs(T_c-T) < eps(typeof(T))
        return (p_c,V_c,V_c)
    end
        #@error "initial temperature $T greater than critical temperature $T_c. returning NaN"
    if 0.6*T_c < T < T_c
        x0 = x0_sat_pure(model,T,crit)
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
