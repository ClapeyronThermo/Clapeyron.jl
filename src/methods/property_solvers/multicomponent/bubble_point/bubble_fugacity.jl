"""
    FugBubblePressure(kwargs...)

Function to compute [`bubble_pressure`](@ref) via fugacity coefficients. First it uses
successive substitution to update the phase composition and a outer newtown
loop to update the pressure. If no convergence is reached after `itmax_newton`
iterations, the system is solved using a multidimensional non-linear
system of equations.

Inputs:
- `y0 = nothing`: optional, initial guess for the vapor phase composition
- `p0 = nothing`: optional, initial guess for the bubble pressure `[Pa]`
- `vol0 = nothing`: optional, initial guesses for the liquid and vapor phase volumes `[m³]`
- `itmax_newton = 10`: optional, number of iterations to update the pressure using newton's method
- `itmax_ss = 5`: optional, number of iterations to update the liquid phase composition using successive substitution
- `tol_x = 1e-8`: optional, tolerance to stop successive substitution cycle
- `tol_p = 1e-8`: optional, tolerance to stop newton cycle
- `tol_of = 1e-8`: optional, tolerance to check if the objective function is zero.
- `nonvolatiles = nothing`: optional, Vector of strings containing non volatile compounds. those will be set to zero on the vapour phase.
- `second_order`: optional, decide if the algorithm uses second order information when updating the guess estimates. Second order methods are slower, but more reliable.
"""
struct FugBubblePressure{T} <: BubblePointMethod
    vol0::Union{Nothing,Tuple{T,T}}
    p0::Union{Nothing,T}
    y0::Union{Nothing,Vector{T}}
    nonvolatiles::Union{Nothing,Vector{String}}
    f_limit::Float64
    atol::Float64
    rtol::Float64
    max_iters::Int
    itmax_newton::Int
    itmax_ss::Int
    tol_y::Float64
    tol_p::Float64
    tol_of::Float64
    second_order::Bool
end

function Solvers.primalval(method::FugBubblePressure{T}) where T
    if T == Nothing
        return Solvers.primalval_struct(method,T)
    else
        return Solvers.primalval_struct(method,Solvers.primal_eltype(T))
    end
end

function FugBubblePressure(;vol0 = nothing,
                                p0 = nothing,
                                y0 = nothing,
                                nonvolatiles = nothing,
                                f_limit = 0.0,
                                atol = 1e-8,
                                rtol = 1e-12,
                                max_iters = 10^4,
                                itmax_newton = 10,
                                itmax_ss = 5,
                                tol_y = 1e-8,
                                tol_p = 1e-8,
                                tol_of = 1e-8,
                                second_order = true)

    if p0 == y0 == vol0 == nothing
        return FugBubblePressure{Nothing}(vol0,p0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_y,tol_p,tol_of,second_order)
    elseif (p0 == y0 == nothing) && !isnothing(vol0)
        vl,vv = promote(vol0[1],vol0[2])
        return FugBubblePressure{typeof(vl)}(vol0,p0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_y,tol_p,tol_of,second_order)
    elseif (vol0 == y0 == nothing) && !isnothing(p0)
        p0 = float(p0)
        return FugBubblePressure{typeof(p0)}(vol0,p0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_y,tol_p,tol_of,second_order)
    elseif (p0 == vol0 == nothing) && !isnothing(y0)
        T = eltype(y0)
        return FugBubblePressure{T}(vol0,p0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_y,tol_p,tol_of,second_order)
    elseif !isnothing(vol0) && !isnothing(p0) && !isnothing(y0)
        vl,vv,p0,_ = promote(vol0[1],vol0[2],p0,first(y0))
        T = eltype(vl)
        y0 = convert(Vector{T},y0)
        return FugBubblePressure{T}(vol0,p0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_y,tol_p,tol_of,second_order)
    elseif !isnothing(vol0) && !isnothing(y0)
        vl,vv,_ = promote(vol0[1],vol0[2],first(y0))
        T = eltype(vl)
        y0 = convert(Vector{T},y0)
        return FugBubblePressure{T}(vol0,p0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_y,tol_p,tol_of,second_order)
    elseif  !isnothing(p0) && !isnothing(y0)
        p0,_ = promote(p0,first(y0))
        T = eltype(p0)
        y0 = convert(Vector{T},y0)
        return FugBubblePressure{T}(vol0,p0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_y,tol_p,tol_of,second_order)
    else
        throw(error("invalid specification for bubble pressure"))
    end
end

function bubble_pressure_impl(model::RestrictedEquilibriaModel,T,x,method::FugBubblePressure)
    wrapper = PTFlashWrapper(model,NaN,T,x,:vle)
    return bubble_pressure_impl(wrapper,T,x,method)
end

function bubble_pressure_impl(model::CompositeModel,T,x, method::FugBubblePressure)
    return bubble_pressure_impl(model.fluid,T,x,method)
end

function bubble_pressure_impl(model::EoSModel, T, x, method::FugBubblePressure)
    nonvolatiles = method.nonvolatiles
    volatiles = comps_in_equilibria(component_list(model),nonvolatiles)
    _vol0,_p0,_y0 = method.vol0,method.p0,method.y0
    p0,volx,voly,y0 = bubble_pressure_init(model,T,x,_vol0,_p0,_y0,volatiles)
    vol0 = (volx,voly)
    model_y,_ = index_reduction(model,volatiles)
    y0 = y0[volatiles]
    y0 ./= sum(y0)

    data = FugData(FugEnum.BUBBLE_PRESSURE,
                    method.itmax_ss,
                    method.itmax_newton,
                    method.tol_p,
                    method.tol_y,
                    method.tol_of,
                    method.second_order,
                    false)

    cache = fug_bubbledew_cache(model,model_y,p0,T,x,y0,Val{false}())

    if all(volatiles)
        converged,res = _fug_OF_ss(model,p0,T,x,y0,vol0,data,cache)
    else
        converged,res = _fug_OF_ss(model,model_y,p0,T,x,y0,vol0,volatiles,data,cache)
    end

    p,T,x,y,vol = res
    volx,voly = vol

    if converged || isnan(volx) || isnan(voly)
        if iszero(volx) && model isa PTFlashWrapper
            vx = volume(model,p,T,x,phase = :l)
            volx = oftype(volx,vx)
        end
        return p,volx,voly,index_expansion(y,volatiles)
    end

    lnK,K,w,w_old,w_calc,w_restart,vol_cache,Hϕx,Hϕy = cache
    inc0 = vcat(lnK, log(p))
    vol_cache[] = (volx,voly)
    opts = NLSolvers.NEqOptions(method)
    if all(volatiles)
        problem = _fug_OF_neq(model,T,x,data,cache)
        sol = Solvers.nlsolve(problem, inc0, Solvers.LineSearch(Solvers.Newton2(inc0)),opts)
        inc = Solvers.x_sol(sol)
        !__check_convergence(sol) && (inc.= NaN)
    else
        problem = _fug_OF_neq(model,model_y,T,x,volatiles,data,cache)
        sol = Solvers.nlsolve(problem, inc0, Solvers.LineSearch(Solvers.Newton2(inc0)),opts)
        inc = Solvers.x_sol(sol)
        !__check_convergence(sol) && (inc.= NaN)
    end

    lnp = inc[end]
    lnK = @view inc[1:(end-1)]
    y_r = exp.(lnK) .* @view(x[volatiles])
    y = index_expansion(y_r,volatiles)
    p = exp(lnp)
    volx,voly = vol_cache[]
    if iszero(volx) && model isa PTFlashWrapper
        vx = volume(model,p,T,x,phase = :l)
        volx = oftype(volx,vx)
    end
    return p, volx, voly, y
end

"""
    FugBubbleTemperature(kwargs...)

Function to compute bubble pressure via fugacity coefficients.
First it uses successive substitution to update the phase composition and a outer newton (or secant) loop to update the temperature.
If no convergence is reached after `itmax_newton` iterations, the system is solved using a multidimensional non-linear system of equations.


Inputs:
- `y = nothing`: optional, initial guess for the vapor phase composition.
- `T0 = nothing`: optional, initial guess for the bubble temperature `[K]`.
- `vol0 = nothing`: optional, initial guesses for the liquid and vapor phase volumes `[m³]`
- `itmax_newton = 10`: optional, number of iterations to update the pressure using newton's (or secant) method
- `itmax_ss = 5`: optional, number of iterations to update the liquid phase composition using successive substitution
- `tol_x = 1e-8`: optional, tolerance to stop successive substitution cycle
- `tol_T = 1e-8`: optional, tolerance to stop newton cycle
- `tol_of = 1e-8`: optional, tolerance to check if the objective function is zero.
- `nonvolatiles`: optional, Vector of strings containing non volatile compounds. those will be set to zero on the vapour phase.
- `second_order`: optional, decide if the algorithm uses second order information when updating the guess estimates. Second order methods are slower, but more reliable.
"""
struct FugBubbleTemperature{T} <: BubblePointMethod
    vol0::Union{Nothing,Tuple{T,T}}
    T0::Union{Nothing,T}
    y0::Union{Nothing,Vector{T}}
    nonvolatiles::Union{Nothing,Vector{String}}
    f_limit::Float64
    atol::Float64
    rtol::Float64
    max_iters::Int
    itmax_newton::Int
    itmax_ss::Int
    tol_y::Float64
    tol_T::Float64
    tol_of::Float64
    second_order::Bool
end

function Solvers.primalval(method::FugBubbleTemperature{T}) where T
    if T == Nothing
        return Solvers.primalval_struct(method,T)
    else
        return Solvers.primalval_struct(method,Solvers.primal_eltype(T))
    end
end

function FugBubbleTemperature(;vol0 = nothing,
    T0 = nothing,
    y0 = nothing,
    nonvolatiles = nothing,
    f_limit = 0.0,
    atol = 1e-8,
    rtol = 1e-12,
    max_iters = 10^4,
    itmax_newton = 10,
    itmax_ss = 5,
    tol_y = 1e-8,
    tol_T = 1e-8,
    tol_of = 1e-8,
    second_order = true)

    if T0 == y0 == vol0 == nothing
        return FugBubbleTemperature{Nothing}(vol0,T0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_y,tol_T,tol_of,second_order)
    elseif (T0 == y0 == nothing) && !isnothing(vol0)
        vl,vv = promote(vol0[1],vol0[2])
        return FugBubbleTemperature{typeof(vl)}(vol0,T0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_y,tol_T,tol_of,second_order)
    elseif (vol0 == y0 == nothing) && !isnothing(T0)
        T0 = float(T0)
        return FugBubbleTemperature{typeof(T0)}(vol0,T0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_y,tol_T,tol_of,second_order)
    elseif (T0 == vol0 == nothing) && !isnothing(y0)
        T = eltype(y0)
        return FugBubbleTemperature{T}(vol0,T0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_y,tol_T,tol_of,second_order)
    elseif !isnothing(vol0) && !isnothing(T0) && !isnothing(y0)
        vl,vv,T0,_ = promote(vol0[1],vol0[2],T0,first(y0))
        T = eltype(vl)
        y0 = convert(Vector{T},y0)
        return FugBubbleTemperature{T}(vol0,T0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_y,tol_T,tol_of,second_order)
    elseif !isnothing(vol0) && !isnothing(y0)
        vl,vv,_ = promote(vol0[1],vol0[2],first(y0))
        T = eltype(vl)
        y0 = convert(Vector{T},y0)
        return FugBubbleTemperature{T}(vol0,T0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_y,tol_T,tol_of,second_order)
    elseif  !isnothing(T0) && !isnothing(y0)
        T0,_ = promote(T0,first(y0))
        T = eltype(T0)
        y0 = convert(Vector{T},y0)
        return FugBubbleTemperature{T}(vol0,T0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_y,tol_T,tol_of,second_order)
    else
        throw(error("invalid specification for bubble temperature"))
    end
end

function bubble_temperature_impl(model::RestrictedEquilibriaModel,p,x,method::FugBubbleTemperature)
    wrapper = PTFlashWrapper(model,p,NaN,x,:vle)
    return bubble_temperature_impl(wrapper,p,x,method)
end

function bubble_temperature_impl(model::CompositeModel,p,x,method::FugBubbleTemperature)
    return bubble_temperature_impl(model.fluid,p,x,method)
end

function bubble_temperature_impl(model::EoSModel, p, x, method::FugBubbleTemperature)
    nonvolatiles = method.nonvolatiles
    volatiles = comps_in_equilibria(component_list(model),nonvolatiles)
    _vol0,_T0,_y0 = method.vol0,method.T0,method.y0
    T0,volx,voly,y0 = bubble_temperature_init(model,p,x,_vol0,_T0,_y0,volatiles)
    vol0 = (volx,voly)
    model_y,_ = index_reduction(model,volatiles)
    y0 = y0[volatiles]
    y0 ./= sum(y0)

    data = FugData(FugEnum.BUBBLE_TEMPERATURE,
                    method.itmax_ss,
                    method.itmax_newton,
                    method.tol_T,
                    method.tol_y,
                    method.tol_of,
                    method.second_order,
                    false)

    cache = fug_bubbledew_cache(model,model_y,p,T0,x,y0,Val{true}())

    if all(volatiles)
        converged,res = _fug_OF_ss(model,p,T0,x,y0,vol0,data,cache)
    else
        converged,res = _fug_OF_ss(model,model_y,p,T0,x,y0,vol0,volatiles,data,cache)
    end

    p,T,x,y,vol = res
    volx,voly = vol
    opts = NLSolvers.NEqOptions(method)
    if converged || isnan(volx) || isnan(voly)
        if iszero(volx) && model isa PTFlashWrapper
            vx = volume(model,p,T,x,phase = :l)
            volx = oftype(volx,vx)
        end
        return T,volx,voly,index_expansion(y,volatiles)
    end

    lnK,K,w,w_old,w_calc,w_restart,vol_cache,Hϕx,Hϕy = cache
    inc0 = vcat(lnK, log(T))
    vol_cache[] = (volx,voly)
    
    if all(volatiles)
        problem = _fug_OF_neq(model,p,x,data,cache)
        sol = Solvers.nlsolve(problem, inc0, Solvers.LineSearch(Solvers.Newton2(inc0)),opts)
        inc = Solvers.x_sol(sol)
        !__check_convergence(sol) && (inc.= NaN)
    else
        problem = _fug_OF_neq(model,model_y,p,x,volatiles,data,cache)
        sol = Solvers.nlsolve(problem, inc0, Solvers.LineSearch(Solvers.Newton2(inc0)),opts)
        inc = Solvers.x_sol(sol)
        !__check_convergence(sol) && (inc.= NaN)
    end

    lnT = inc[end]
    lnK = @view inc[1:(end-1)]
    y_r = exp.(lnK) .* @view(x[volatiles])
    y = index_expansion(y_r,volatiles)
    T = exp(lnT)
    volx,voly = vol_cache[]
    if iszero(volx) && model isa PTFlashWrapper
        vx = volume(model,p,T,x,phase = :l)
        volx = oftype(volx,vx)
    end
    return T, volx, voly, y
end

export FugBubblePressure, FugBubbleTemperature