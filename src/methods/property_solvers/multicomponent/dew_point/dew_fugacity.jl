########## Dew pressure calculation
"""
    OF_dewPx!(model::EoSModel,modely y, T, vol_cache)
    OF_dewPx!(modelx::EoSModel,modely::EoSModel y, T, vol_cache,_views)

Objective function to compute dew pressure using a multidimensional
system of equations via fugacity coefficients.

Inputs:
- `model`: general equation of state model
- `modelx`: liquid equation of state model, if any noncondensable compounds are present
- `modely`: vapour equation of state model
- `y`: vapour phase composition
- `T`: temperature `[K]`
- `vol_cache`: array used to update the phases' volumes
_ `condensable`: condensable component indices, if any noncondensable compounds are present

Returns: NLSolvers.NEqProblem
"""
function OF_dewPx! end

function OF_dewPx!(model, y, T, vol_cache)
    return _fug_OF_neqsystem(model, nothing, y, nothing, T, vol_cache, FugEnum.DEW_PRESSURE, (:liquid,:vapor))
end

function OF_dewPx!(modelx,modely, y, T, vol_cache,condensable)
    return _fug_OF_neqsystem(modelx, modely, nothing, y, nothing, T, vol_cache, FugEnum.DEW_PRESSURE, (:liquid,:vapor), condensable)
end

"""
    function dew_pressure_fug(model::EoSModel, T, y, x0, p0; vol0=(nothing,nothing),
                             itmax_newton = 10, itmax_ss = 5, tol_x = 1e-8,
                             tol_p = 1e-8, tol_of = 1e-8)

Function to compute dew pressure via fugacity coefficients. First it uses
successive substitution to update the phase composition and a outer newtown
loop to update the pressure. If no convergence is reached after itmax_newton
iterations, the system is solved using a multidimensional non-linear
systems of equations.

Inputs:
- `model`: equation of state model
- `T`: dew temperature `[K]`
- `y`: vapor phase composition
- `x0`: initial guess for the liquid phase composition `[m³]`
- `p0`: initial guess for the dew pressure `[Pa]`
- `vol0`: optional, initial guesses for the liquid and vapor phase volumes
- `itmax_newton`: optional, number of iterations to update the pressure using newton's method
- `itmax_ss`: optional, number of iterations to update the liquid phase composition using successive substitution
- `tol_x`: optional, tolerance to stop successive substitution cycle
- `tol_p`: optional, tolerance to stop newton cycle
- `tol_of`: optional, tolerance to check if the objective function is zero.
- `noncondensables`: optional, Vector of strings containing non condensable compounds. those will be set to zero on the liquid phase.

Returns:
- `p`: dew pressure `[Pa]`
- `volx`: saturared liquid volume `[m³]`
- `voly`: saturared vapor volume `[m³]`
- `x`: saturated liquid composition
"""
function dew_pressure_fug(model::EoSModel, T, y, x0, p0; vol0=(nothing,nothing),
                             itmax_newton = 10, itmax_ss = 5, tol_x = 1e-8,
                             tol_p = 1e-8, tol_of = 1e-8, noncondensables = nothing)

     # Setting the initial guesses for volumes
    vol0 === nothing && (vol0 = (nothing,nothing))
    volx, voly = vol0

    #check if noncondensables are set
    if !isnothing(noncondensables)
        condensables = [!in(x,noncondensables) for x in model.components]
        model_x,condensables = index_reduction(model,condensables)
        x0 = x0[condensables]
        x0 = x0/sum(x0)
    else
        condensables = fill(true,length(model))
        model_x = nothing
    end

    converged,res = _fug_OF_ss(model_x,model,p0,T,x0,y,vol0,FugEnum.DEW_PRESSURE,condensables;itmax_ss = itmax_ss, itmax_newton = itmax_newton,tol_pT = tol_p, tol_xy = tol_x, tol_of = tol_of)
    p,T,x,y,vol,lnK = res
    volx,voly = vol
    if converged
        return p,volx,voly,index_expansion(x,condensables)
    else
        inc0 = vcat(lnK, log(p))
        vol_cache = [volx, voly]
        problem = OF_dewPx!(model_x,model, y, T, vol_cache,condensables)
        sol = Solvers.nlsolve(problem, inc0, Solvers.LineSearch(Solvers.Newton2(inc0)))
        inc = Solvers.x_sol(sol)
        !all(<(sol.options.f_abstol),sol.info.best_residual) && (inc .= NaN)
        lnK = inc[1:(end-1)]
        lnp = inc[end]

        x_r = y[condensables] ./ exp.(lnK)
        x = index_expansion(x_r,condensables)
        p = exp(lnp)
        volx, voly = vol_cache
     end

     return p, volx, voly, x
end
"""
    FugDewPressure(kwargs...)

Method to compute [`dew_pressure`](@ref) via fugacity coefficients. First it uses
successive substitution to update the phase composition and a outer newtown
loop to update the pressure. If no convergence is reached after `itmax_newton`
iterations, the system is solved using a multidimensional non-linear
system of equations.

Inputs:
- `x0 = nothing`: optional, initial guess for the liquid phase composition
- `p0 = nothing`: optional, initial guess for the dew pressure `[Pa]`
- `vol0 = nothing`: optional, initial guesses for the liquid and vapor phase volumes `[m³]`
- `itmax_newton = 10`: optional, number of iterations to update the pressure using newton's method
- `itmax_ss = 5`: optional, number of iterations to update the liquid phase composition using successive substitution
- `tol_x = 1e-8`: optional, tolerance to stop successive substitution cycle
- `tol_p = 1e-8`: optional, tolerance to stop newton cycle
- `tol_of = 1e-8`: optional, tolerance to check if the objective function is zero.
- `noncondensables = nothing`: optional, Vector of strings containing non condensable compounds. those will be set to zero on the liquid phase.
"""
struct FugDewPressure{T} <: DewPointMethod
    vol0::Union{Nothing,Tuple{T,T}}
    p0::Union{Nothing,T}
    x0::Union{Nothing,Vector{T}}
    noncondensables::Union{Nothing,Vector{String}}
    f_limit::Float64
    atol::Float64
    rtol::Float64
    max_iters::Int
    itmax_newton::Int
    itmax_ss::Int
    tol_x::Float64
    tol_p::Float64
    tol_of::Float64
end


function FugDewPressure(;vol0 = nothing,
                                p0 = nothing,
                                x0 = nothing,
                                noncondensables = nothing,
                                f_limit = 0.0,
                                atol = 1e-8,
                                rtol = 1e-12,
                                max_iters = 10^4,
                                itmax_newton = 10,
                                itmax_ss = 5,
                                tol_x = 1e-8,
                                tol_p = 1e-8,
                                tol_of = 1e-8)


    if p0 == x0 == vol0 == nothing
        return FugDewPressure{Nothing}(vol0,p0,x0,noncondensables,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_x,tol_p,tol_of)
    elseif (p0 == x0 == nothing) && !isnothing(vol0)
        vl,vv = promote(vol0[1],vol0[2])
        return FugDewPressure{typeof(vl)}(vol0,p0,x0,noncondensables,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_x,tol_p,tol_of)
    elseif (vol0 == x0 == nothing) && !isnothing(p0)
        p0 = float(p0)
        return FugDewPressure{typeof(p0)}(vol0,p0,x0,noncondensables,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_x,tol_p,tol_of)
    elseif (p0 == vol0 == nothing) && !isnothing(x0)
        T = eltype(x0)
        return FugDewPressure{T}(vol0,p0,x0,noncondensables,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_x,tol_p,tol_of)
    elseif !isnothing(vol0) && !isnothing(p0) && !isnothing(x0)
        vl,vv,p0,_ = promote(vol0[1],vol0[2],p0,first(x0))
        T = eltype(vl)
        x0 = convert(Vector{T},x0)
        return FugDewPressure{T}(vol0,p0,x0,noncondensables,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_x,tol_p,tol_of)
    elseif !isnothing(vol0) && !isnothing(x0)
        vl,vv,_ = promote(vol0[1],vol0[2],first(x0))
        T = eltype(vl)
        x0 = convert(Vector{T},x0)
        return FugDewPressure{T}(vol0,p0,x0,noncondensables,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_x,tol_p,tol_of)
    elseif  !isnothing(p0) && !isnothing(x0)
        p0,_ = promote(p0,first(x0))
        T = eltype(p0)
        x0 = convert(Vector{T},x0)
        return FugDewPressure{T}(vol0,p0,x0,noncondensables,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_x,tol_p,tol_of)
    else
        throw(error("invalid specification for dew pressure"))
    end
end


function dew_pressure_impl(model::EoSModel, T, y ,method::FugDewPressure)
    
    if !isnothing(method.noncondensables)
        condensables = [!in(x,method.noncondensables) for x in model.components]
    else
        condensables = fill(true,length(model))
    end

    _vol0,_p0,_x0 = method.vol0,method.p0,method.x0
    p0,vl,vv,x0 = dew_pressure_init(model,T,y,_vol0,_p0,_x0,condensables)
    itmax_newton = method.itmax_newton
    itmax_ss = method.itmax_ss
    tol_x = method.tol_x
    tol_p = method.tol_p
    tol_of = method.tol_of
    vol0 = (vl,vv)
    noncondensables = method.noncondensables
    return dew_pressure_fug(model,T,y,x0,p0;vol0,itmax_newton,itmax_ss,tol_x,tol_p,tol_of,noncondensables)
end

################# Dew temperature calculation

"""
    OF_dewTx!(model::EoSModel, y, p, vol_cache)
    OF_dewTx!(modelx::EoSModel,modely::EoSModel y, p, vol_cache,_views)

Objective function to compute dew temperature using a multidimensional
system of equations via fugacity coefficients.

Inputs:
- `model`: general equation of state model
- `modelx`: liquid equation of state model, if any noncondensable compounds are present
- `modely`: vapour equation of state model
- `P`: pressure `[Pa]`
- `vol_cache`: array used to update the phases' volumes
_ `condensable`: condensable component indices, if any noncondensable compounds are present

Returns: NLSolvers.NEqProblem
"""
function OF_dewTx! end

function OF_dewTx!(model, y, p, vol_cache)
    return _fug_OF_neqsystem(model, nothing, y, p, nothing, vol_cache, FugEnum.DEW_TEMPERATURE, (:liquid,:vapor))
end

function OF_dewTx!(model,modely, y, p, vol_cache,condensable)
    return _fug_OF_neqsystem(model, modely, nothing, y, p, nothing, vol_cache, FugEnum.DEW_TEMPERATURE, (:liquid,:vapor), condensable)
end

"""
    dew_temperature_fug(model::EoSModel, p, y, x0, T0; vol0=(nothing,nothing),
                             itmax_newton = 10, itmax_ss = 5, tol_x = 1e-8,
                             tol_T = 1e-8, tol_of = 1e-8)

Function to compute dew temperature via fugacity coefficients. First it uses
successive substitution to update the phase composition and a outer newtown
loop to update the temperature. If no convergence is reached after
itmax_newton iterations, the system is solved using a multidimensional
non-linear system of equations.

Inputs:
model: equation of state model
- `P`: pressure `[Pa]`
- `y`: vapor phase composition
- `x0`: initial guess for the liquid phase composition
- `T0`: initial guess for the dew temperature `[K]`
- `vol0`: optional, initial guesses for the liquid and vapor phase volumes `[m³]`
- `itmax_newton`: optional, number of iterations to update the temperature using newton's method
- `itmax_ss`: optional, number of iterations to update the liquid phase composition using successive substitution
- `tol_x`: optional, tolerance to stop successive substitution cycle
- `tol_T`: optional, tolerance to stop newton cycle
- `tol_of`: optional, tolerance to check if the objective function is zero.
- `noncondensables`: optional, Vector of strings containing non condensable compounds. those will be set to zero on the liquid phase.

Returns:
`T`: dew temperature `[K]`
`volx`: saturared liquid volume `[m³]`
`voly`: saturared vapor volume `[m³]`
`x`: saturated liquid composition
"""
function dew_temperature_fug(model::EoSModel, p, y, x0, T0; vol0=(nothing,nothing),
                             itmax_newton = 10, itmax_ss = 5, tol_x = 1e-8,
                             tol_T = 1e-8, tol_of = 1e-8,noncondensables = nothing)

     # Setting the initial guesses for volumes
    vol0 === nothing && (vol0 = (nothing,nothing))
    volx, voly = vol0
    #check if noncondensables are set
    if !isnothing(noncondensables)
        condensables = [!in(x,noncondensables) for x in model.components]
        model_x,condensables = index_reduction(model,condensables)
        x0 = x0[condensables]
        x0 = x0/sum(x0)
    else
        condensables = fill(true,length(model))
        model_x = nothing
    end

    converged,res = _fug_OF_ss(model_x,model,p,T0,x0,y,vol0,FugEnum.DEW_TEMPERATURE,condensables;itmax_ss = itmax_ss, itmax_newton = itmax_newton, tol_pT = tol_T, tol_xy = tol_x, tol_of = tol_of)
    p,T,x,y,vol,lnK = res
    volx,voly = vol
    if converged
        return T,volx,voly,index_expansion(x,condensables)
    else
        inc0 = vcat(lnK, log(T))
        vol_cache = [volx, voly]
        problem = OF_dewTx!(model_x,model, y, p, vol_cache,condensables)
        sol = Solvers.nlsolve(problem, inc0, Solvers.LineSearch(Solvers.Newton2(inc0)))
        inc = Solvers.x_sol(sol)
        !all(<(sol.options.f_abstol),sol.info.best_residual) && (inc .= NaN)
        lnK = inc[1:(end-1)]
        lnT = inc[end]
        x_r = y[condensables]./ exp.(lnK)
        x = index_expansion(x_r,condensables)
        T = exp(lnT)
        volx, voly = vol_cache[:]
    end
    return T, volx, voly, x
end

"""
    FugDewTemperature(kwargs...)

Method to compute [`dew_temperature`](@ref) via fugacity coefficients. First it uses
successive substitution to update the phase composition and a outer newtown
loop to update the temperature. If no convergence is reached after
`itmax_newton` iterations, the system is solved using a multidimensional
non-linear system of equations.

Inputs:
- `x0 = nothing`: optional, initial guess for the liquid phase composition
- `T0 = nothing`: optional, initial guess for the dew temperature `[K]`
- `vol0 = nothing`: optional, initial guesses for the liquid and vapor phase volumes `[m³]`
- `itmax_newton = 10`: optional, number of iterations to update the temperature using newton's method
- `itmax_ss = 5`: optional, number of iterations to update the liquid phase composition using successive substitution
- `tol_x = 1e-8`: optional, tolerance to stop successive substitution cycle
- `tol_T = 1e-8`: optional, tolerance to stop newton cycle
- `tol_of = 1e-8`: optional, tolerance to check if the objective function is zero.
- `noncondensables = nothing`: optional, Vector of strings containing non condensable compounds. those will be set to zero on the liquid phase.

"""
struct FugDewTemperature{T} <: DewPointMethod
    vol0::Union{Nothing,Tuple{T,T}}
    T0::Union{Nothing,T}
    x0::Union{Nothing,Vector{T}}
    noncondensables::Union{Nothing,Vector{String}}
    f_limit::Float64
    atol::Float64
    rtol::Float64
    max_iters::Int
    itmax_newton::Int
    itmax_ss::Int
    tol_x::Float64
    tol_T::Float64
    tol_of::Float64
end

function FugDewTemperature(;vol0 = nothing,
    T0 = nothing,
    x0 = nothing,
    noncondensables = nothing,
    f_limit = 0.0,
    atol = 1e-8,
    rtol = 1e-12,
    max_iters = 10^4,
    itmax_newton = 10,
    itmax_ss = 5,
    tol_x = 1e-8,
    tol_T = 1e-8,
    tol_of = 1e-8)

    if T0 == x0 == vol0 == nothing
        return FugDewTemperature{Nothing}(vol0,T0,x0,noncondensables,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_x,tol_T,tol_of)
    elseif (T0 == x0 == nothing) && !isnothing(vol0)
        vl,vv = promote(vol0[1],vol0[2])
        return FugDewTemperature{typeof(vl)}(vol0,T0,x0,noncondensables,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_x,tol_T,tol_of)
    elseif (vol0 == x0 == nothing) && !isnothing(T0)
        T0 = float(T0)
        return FugDewTemperature{typeof(T0)}(vol0,T0,x0,noncondensables,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_x,tol_T,tol_of)
    elseif (T0 == vol0 == nothing) && !isnothing(x0)
        T = eltype(x0)
        return FugDewTemperature{T}(vol0,T0,x0,noncondensables,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_x,tol_T,tol_of)
    elseif !isnothing(vol0) && !isnothing(T0) && !isnothing(x0)
        vl,vv,T0,_ = promote(vol0[1],vol0[2],T0,first(x0))
        T = eltype(vl)
        x0 = convert(Vector{T},x0)
        return FugDewTemperature{T}(vol0,T0,x0,noncondensables,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_x,tol_T,tol_of)
    elseif !isnothing(vol0) && !isnothing(x0)
        vl,vv,_ = promote(vol0[1],vol0[2],first(x0))
        T = eltype(vl)
        x0 = convert(Vector{T},x0)
        return FugDewTemperature{T}(vol0,T0,x0,noncondensables,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_x,tol_T,tol_of)
    elseif  !isnothing(T0) && !isnothing(x0)
        T0,_ = promote(T0,first(x0))
        T = eltype(T0)
        x0 = convert(Vector{T},x0)
        return FugDewTemperature{T}(vol0,T0,x0,noncondensables,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_x,tol_T,tol_of)
    else
        throw(error("invalid specification for bubble temperature"))
    end
end

function dew_temperature_impl(model::EoSModel, p, y, method::FugDewTemperature)
    if !isnothing(method.noncondensables)
        condensables = [!in(x,method.noncondensables) for x in model.components]
    else
        condensables = fill(true,length(model))
    end

    _vol0,_T0,_x0 = method.vol0,method.T0,method.x0
    T0,vl,vv,x0 = dew_temperature_init(model,p,y,_vol0,_T0,_x0,condensables)
    itmax_newton = method.itmax_newton
    itmax_ss = method.itmax_ss
    tol_x = method.tol_x
    tol_T = method.tol_T
    tol_of = method.tol_of
    vol0 = (vl,vv)

    noncondensables = method.noncondensables
    return dew_temperature_fug(model,p,y,x0,T0;vol0,itmax_newton,itmax_ss,tol_x,tol_T,tol_of,noncondensables)
end

export FugDewPressure, FugDewTemperature