"""
    OF_bubblepy!(model::EoSModel,modely x, T, vol_cache)
    OF_bubblepy!(modelx::EoSModel,modely::EoSModel x, T, vol_cache,volatile)

Objective function to compute bubble pressure using a multidimensional
system of equations via fugacity coefficients.

Inputs:
- `model`: general equation of state model
- `modelx`: liquid equation of state model
- `modely`: vapour equation of state model, if any nonvolatile compounds are present
- `x`: liquid phase composition
- `T`: temperature [`K`]
- `vol_cache`: array used to update the phases' volumes
_ `volatile`: volatile component indices, if any nonvolatile compounds are present

Returns: NLSolvers.NEqProblem
"""
function OF_bubblepy! end

function OF_bubblepy!(model, x, T, vol_cache)
    return _fug_OF_neqsystem(model,x, nothing, nothing, T, vol_cache,FugEnum.BUBBLE_PRESSURE,(:liquid,:vapor))
end

function OF_bubblepy!(model,modely, x, T, vol_cache,volatile)
    return _fug_OF_neqsystem(model, modely, x, nothing, nothing, T, vol_cache, FugEnum.BUBBLE_PRESSURE, (:liquid,:vapor), volatile)
end

"""
    bubble_pressure_fug(model::EoSModel, T, x, y0, p0; vol0=(nothing,nothing),
                         itmax_newton = 10, itmax_ss = 5, tol_y = 1e-8,
                         tol_p = 1e-8, tol_of = 1e-8)

Function to compute bubble pressure via fugacity coefficients. First it uses
successive substitution to update the phase composition and a outer newtown
loop to update the pressure. If no convergence is reached after itmax_newton
iterations, the system is solved using a multidimensional non-linear
systems of equations.

Inputs:
model: equation of state model
- `T`: bubble temperature [`K`]
- `x`: liquid phase composition
- `y0`: initial guess for the vapor phase composition
- `p0`: initial guess for the bubble pressure [`Pa`]
- `vol0`: optional, initial guesses for the liquid and vapor phase volumes
- `itmax_newton`: optional, number of iterations to update the pressure using newton's method
- `itmax_ss`: optional, number of iterations to update the liquid phase composition using successive substitution
- `tol_x`: optional, tolerance to stop successive substitution cycle
- `tol_p`: optional, tolerance to stop newton cycle
- `tol_of`: optional, tolerance to check if the objective function is zero.
- `nonvolatiles = nothing`: optional, Vector of strings containing non volatile compounds. those will be set to zero on the vapour phase.

Returns:
`p`: bubble pressure
`volx`: saturared liquid volume
`voly`: saturared vapor volume
`y`: saturated vapor composition
"""
function bubble_pressure_fug(model::EoSModel, T, x, y0, p0; vol0=(nothing,nothing),
                         itmax_newton = 10, itmax_ss = 5, tol_y = 1e-8,
                         tol_p = 1e-8, tol_of = 1e-8,nonvolatiles = nothing)


    # Setting the initial guesses for volumes
    vol0 === nothing && (vol0 = (nothing,nothing))
    volx, voly = vol0

    #check if nonvolatiles are set
    if !isnothing(nonvolatiles)
        volatiles = [!in(x,nonvolatiles) for x in model.components]
        model_y,volatiles = index_reduction(model,volatiles)
        y0 = y0[volatiles]
        y0 = y0/sum(y0)
    else
        volatiles = fill(true,length(model))
        model_y = nothing
    end

    converged,res = _fug_OF_ss(model,model_y,p0,T,x,y0,vol0,FugEnum.BUBBLE_PRESSURE,volatiles;itmax_ss = itmax_ss, itmax_newton = itmax_newton,tol_pT = tol_p,tol_xy = tol_y,tol_of=tol_of)
    p,T,x,y,vol,lnK = res
    volx,voly = vol
    if converged
        return p,volx,voly,index_expansion(y,volatiles)
    elseif isnan(volx) || isnan(voly)
        return p,volx,voly,index_expansion(y,volatiles)
    else
        inc0 = vcat(lnK, log(p))
        vol_cache = [volx, voly]
        problem = OF_bubblepy!(model,model_y, x, T, vol_cache,volatiles)
        sol = Solvers.nlsolve(problem, inc0, Solvers.LineSearch(Solvers.Newton2(inc0)))
        inc = Solvers.x_sol(sol)
        !all(<(sol.options.f_abstol),sol.info.best_residual) && (inc .= NaN)
        lnp = inc[end]
        lnK = inc[1:(end-1)]

        y_r = exp.(lnK) .* x[volatiles]
        y = index_expansion(y_r,volatiles)
        p = exp.(lnp)
        volx, voly = vol_cache[:]
        # println("Second order method ", p, " ", y)
    end
    #@show p,volx,voly,y
    return p, volx, voly, y
end

"""
    FugBubblePressure(kwargs...)

Function to compute [`bubble_pressure`](@ref) via fugacity coefficients. First it uses
successive substitution to update the phase composition and a outer newtown
loop to update the pressure. If no convergence is reached after `itmax_newton`
iterations, the system is solved using a multidimensional non-linear
system of equations.

Inputs:
- `y0 = nothing`: optional, initial guess for the vapor phase composition
- `p0 = nothing`: optional, initial guess for the bubble pressure [`Pa`]
- `vol0 = nothing`: optional, initial guesses for the liquid and vapor phase volumes
- `itmax_newton = 10`: optional, number of iterations to update the pressure using newton's method
- `itmax_ss = 5`: optional, number of iterations to update the liquid phase composition using successive substitution
- `tol_x = 1e-8`: optional, tolerance to stop successive substitution cycle
- `tol_p = 1e-8`: optional, tolerance to stop newton cycle
- `tol_of = 1e-8`: optional, tolerance to check if the objective function is zero.
- `nonvolatiles = nothing`: optional, Vector of strings containing non volatile compounds. those will be set to zero on the vapour phase.
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
                                tol_of = 1e-8)

    if p0 == y0 == vol0 == nothing
        return FugBubblePressure{Nothing}(vol0,p0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_y,tol_p,tol_of)
    elseif (p0 == y0 == nothing) && !isnothing(vol0)
        vl,vv = promote(vol0[1],vol0[2])
        return FugBubblePressure{typeof(vl)}(vol0,p0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_y,tol_p,tol_of)
    elseif (vol0 == y0 == nothing) && !isnothing(p0)
        p0 = float(p0)
        return FugBubblePressure{typeof(p0)}(vol0,p0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_y,tol_p,tol_of)
    elseif (p0 == vol0 == nothing) && !isnothing(y0)
        T = eltype(y0)
        return FugBubblePressure{T}(vol0,p0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_y,tol_p,tol_of)
    elseif !isnothing(vol0) && !isnothing(p0) && !isnothing(y0)
        vl,vv,p0,_ = promote(vol0[1],vol0[2],p0,first(y0))
        T = eltype(vl)
        y0 = convert(Vector{T},y0)
        return FugBubblePressure{T}(vol0,p0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_y,tol_p,tol_of)
    elseif !isnothing(vol0) && !isnothing(y0)
        vl,vv,_ = promote(vol0[1],vol0[2],first(y0))
        T = eltype(vl)
        y0 = convert(Vector{T},y0)
        return FugBubblePressure{T}(vol0,p0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_y,tol_p,tol_of)
    elseif  !isnothing(p0) && !isnothing(y0)
        p0,_ = promote(p0,first(y0))
        T = eltype(p0)
        y0 = convert(Vector{T},y0)
        return FugBubblePressure{T}(vol0,p0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_y,tol_p,tol_of)
    else
        throw(error("invalid specification for bubble pressure"))
    end
end

function bubble_pressure_impl(model::EoSModel, T, x,method::FugBubblePressure)
    if !isnothing(method.nonvolatiles)
        volatiles = [!in(x,method.nonvolatiles) for x in model.components]
    else
        volatiles = fill(true,length(model))
    end
    _vol0,_p0,_y0 = method.vol0,method.p0,method.y0
    p0,vl,vv,y0 = bubble_pressure_init(model,T,x,_vol0,_p0,_y0,volatiles)
    itmax_newton = method.itmax_newton
    itmax_ss = method.itmax_ss
    tol_y = method.tol_y
    tol_p = method.tol_p
    tol_of = method.tol_of
    vol0 = (vl,vv)
    nonvolatiles = method.nonvolatiles
    return bubble_pressure_fug(model,T,x,y0,p0;vol0,itmax_newton,itmax_ss,tol_y,tol_p,tol_of,nonvolatiles)
end

################ Bubble temperature solver

"""
    OF_bubbleTy!(model::EoSModel, y, p, vol_cache)
    OF_bubblepy!(modelx::EoSModel,modely::EoSModel x, p, vol_cache,_views)

Objective function to compute bubble temperature using a multidimensional
system of equations via fugacity coefficients.

Inputs:
- `model`: general equation of state model
- `modelx`: liquid equation of state model
- `modely`: vapour equation of state model, if any nonvolatile compounds are present
- `x`: liquid phase composition
- `p`: pressure [`Pa`]
- `vol_cache`: array used to update the phases' volumes
_ `volatile`: volatile component indices, if any nonvolatile compounds are present

Returns: NLSolvers.NEqProblem
"""
function OF_bubbleTy! end

function OF_bubbleTy!(model, x, p, vol_cache)
    return _fug_OF_neqsystem(model,x, nothing, p, nothing, vol_cache, FugEnum.BUBBLE_TEMPERATURE, (:liquid,:vapor))
end

function OF_bubbleTy!(model,modely, x, p, vol_cache,volatile)
    return _fug_OF_neqsystem(model, modely, x, nothing, p, nothing, vol_cache, FugEnum.BUBBLE_TEMPERATURE, (:liquid,:vapor), volatile)
end

"""
    bubble_temperature_fug(model::EoSModel, p, x, y0, T0; vol0=(nothing,nothing),
                         itmax_newton = 10, itmax_ss = 5, tol_y = 1e-8,
                         tol_T = 1e-8, tol_of = 1e-8)

Function to compute bubble temperature via fugacity coefficients. First it uses
successive substitution to update the phase composition and a outer newtown
loop to update the temperature. If no convergence is reached after
itmax_newton iterations, the system is solved using a multidimensional
non-linear systems of equations.

Inputs:
- model: equation of state model
- `P`: pressure [`Pa`]
- `x`: liquid phase composition
- `y`: initial guess for the vapor phase composition
- `T0`: initial guess for the bubble temperature [`K`]
- `vol0`: optional, initial guesses for the liquid and vapor phase volumes
- `itmax_newton`: optional, number of iterations to update the temperature using newton's method
- `itmax_ss`: optional, number of iterations to update the liquid phase composition using successive substitution
- `tol_x`: optional, tolerance to stop successive substitution cycle
- `tol_T`: optional, tolerance to stop newton cycle
- `tol_of`: optional, tolerance to check if the objective function is zero.
- `nonvolatiles`: optional, Vector of strings containing non volatile compounds. those will be set to zero on the vapour phase.

Returns:
- `T`: bubble temperature
- `volx`: saturared liquid volume
- `voly`: saturared vapor volume
- `y`: saturated vapor composition
"""
function bubble_temperature_fug(model::EoSModel, p, x, y0, T0; vol0=(nothing,nothing),
                         itmax_newton = 10, itmax_ss = 5, tol_y = 1e-8,
                         tol_T = 1e-8, tol_of = 1e-8,nonvolatiles = nothing)

    # Setting the initial guesses for volumes
    
    vol0 === nothing && (vol0 = (nothing,nothing))
    volx, voly = vol0

    #check if nonvolatiles are set
    if !isnothing(nonvolatiles)
        volatiles = [!in(x,nonvolatiles) for x in model.components]
        model_y,volatiles = index_reduction(model,volatiles)
        y0 = y0[volatiles]
        y0 = y0/sum(y0)
    else
        volatiles = fill(true,length(model))
        model_y = nothing
    end

    converged,res = _fug_OF_ss(model,model_y,p,T0,x,y0,vol0,FugEnum.BUBBLE_TEMPERATURE,volatiles;itmax_ss = itmax_ss, itmax_newton = itmax_newton, tol_pT = tol_T, tol_xy = tol_y, tol_of = tol_of)
    p,T,x,y,vol,lnK = res
    volx,voly = vol
    if converged
        return T,volx,voly,index_expansion(y,volatiles)
    else
        inc0 = vcat(lnK, log(T))
        vol_cache = [volx, voly]
        problem = OF_bubbleTy!(model,model_y, x, p, vol_cache,volatiles)
        sol = Solvers.nlsolve(problem, inc0, Solvers.LineSearch(Solvers.Newton2(inc0)))
        inc = Solvers.x_sol(sol)
        !all(<(sol.options.f_abstol),sol.info.best_residual) && (inc .= NaN)
        lnK = inc[1:(end-1)]
        lnT = inc[end]

        y_r = exp.(lnK) .* x[volatiles]
        y = index_expansion(y_r,volatiles)
        T = exp(lnT)
        volx, voly = vol_cache[:]
    end

    return T, volx, voly, y
end
"""
    FugBubbleTemperature(kwargs...)

Method to compute [`bubble_temperature`](@ref) via fugacity coefficients. First it uses
successive substitution to update the phase composition and a outer newtown
loop to update the temperature. If no convergence is reached after
`itmax_newton` iterations, the system is solved using a multidimensional
non-linear system of equations.

Inputs:
- `y = nothing`: optional, initial guess for the vapor phase composition.
- `T0 = nothing`: optional, initial guess for the bubble temperature [`K`].
- `vol0 = nothing`: optional, initial guesses for the liquid and vapor phase volumes
- `itmax_newton = 10`: optional, number of iterations to update the temperature using newton's method
- `itmax_ss = 5`: optional, number of iterations to update the liquid phase composition using successive substitution
- `tol_x = 1e-8`: optional, tolerance to stop successive substitution cycle
- `tol_T = 1e-8`: optional, tolerance to stop newton cycle
- `tol_of = 1e-8`: optional, tolerance to check if the objective function is zero.
- `nonvolatiles`: optional, Vector of strings containing non volatile compounds. those will be set to zero on the vapour phase.
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
    tol_of = 1e-8)

    if T0 == y0 == vol0 == nothing
        return FugBubbleTemperature{Nothing}(vol0,T0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_y,tol_T,tol_of)
    elseif (T0 == y0 == nothing) && !isnothing(vol0)
        vl,vv = promote(vol0[1],vol0[2])
        return FugBubbleTemperature{typeof(vl)}(vol0,T0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_y,tol_T,tol_of)
    elseif (vol0 == y0 == nothing) && !isnothing(T0)
        T0 = float(T0)
        return FugBubbleTemperature{typeof(T0)}(vol0,T0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_y,tol_T,tol_of)
    elseif (T0 == vol0 == nothing) && !isnothing(y0)
        T = eltype(y0)
        return FugBubbleTemperature{T}(vol0,T0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_y,tol_T,tol_of)
    elseif !isnothing(vol0) && !isnothing(T0) && !isnothing(y0)
        vl,vv,T0,_ = promote(vol0[1],vol0[2],T0,first(y0))
        T = eltype(vl)
        y0 = convert(Vector{T},y0)
        return FugBubbleTemperature{T}(vol0,T0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_y,tol_T,tol_of)
    elseif !isnothing(vol0) && !isnothing(y0)
        vl,vv,_ = promote(vol0[1],vol0[2],first(y0))
        T = eltype(vl)
        y0 = convert(Vector{T},y0)
        return FugBubbleTemperature{T}(vol0,T0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_y,tol_T,tol_of)
    elseif  !isnothing(T0) && !isnothing(y0)
        T0,_ = promote(T0,first(y0))
        T = eltype(T0)
        y0 = convert(Vector{T},y0)
        return FugBubbleTemperature{T}(vol0,T0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_y,tol_T,tol_of)
    else
        throw(error("invalid specification for bubble temperature"))
    end
end

function bubble_temperature_impl(model::EoSModel, p, x, method::FugBubbleTemperature)
    if !isnothing(method.nonvolatiles)
        volatiles = [!in(x,method.nonvolatiles) for x in model.components]
    else
        volatiles = fill(true,length(model))
    end
    _vol0,_T0,_y0 = method.vol0,method.T0,method.y0
    T0,vl,vv,y0 = bubble_temperature_init(model,p,x,_vol0,_T0,_y0,volatiles)
    itmax_newton = method.itmax_newton
    itmax_ss = method.itmax_ss
    tol_y = method.tol_y
    tol_T = method.tol_T
    tol_of = method.tol_of
    vol0 = (vl,vv)
    nonvolatiles = method.nonvolatiles
    return bubble_temperature_fug(model,p,x,y0,T0;vol0,itmax_newton,itmax_ss,tol_y,tol_T,tol_of,nonvolatiles)
end

export FugBubblePressure, FugBubbleTemperature