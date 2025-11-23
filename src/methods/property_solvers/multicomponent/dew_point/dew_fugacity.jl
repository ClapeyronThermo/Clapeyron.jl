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
    second_order::Bool
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
                                tol_of = 1e-8,
                                second_order = true)


    if p0 == x0 == vol0 == nothing
        return FugDewPressure{Nothing}(vol0,p0,x0,noncondensables,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_x,tol_p,tol_of,second_order)
    elseif (p0 == x0 == nothing) && !isnothing(vol0)
        vl,vv = promote(vol0[1],vol0[2])
        return FugDewPressure{typeof(vl)}(vol0,p0,x0,noncondensables,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_x,tol_p,tol_of,second_order)
    elseif (vol0 == x0 == nothing) && !isnothing(p0)
        p0 = float(p0)
        return FugDewPressure{typeof(p0)}(vol0,p0,x0,noncondensables,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_x,tol_p,tol_of,second_order)
    elseif (p0 == vol0 == nothing) && !isnothing(x0)
        T = eltype(x0)
        return FugDewPressure{T}(vol0,p0,x0,noncondensables,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_x,tol_p,tol_of,second_order)
    elseif !isnothing(vol0) && !isnothing(p0) && !isnothing(x0)
        vl,vv,p0,_ = promote(vol0[1],vol0[2],p0,first(x0))
        T = eltype(vl)
        x0 = convert(Vector{T},x0)
        return FugDewPressure{T}(vol0,p0,x0,noncondensables,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_x,tol_p,tol_of,second_order)
    elseif !isnothing(vol0) && !isnothing(x0)
        vl,vv,_ = promote(vol0[1],vol0[2],first(x0))
        T = eltype(vl)
        x0 = convert(Vector{T},x0)
        return FugDewPressure{T}(vol0,p0,x0,noncondensables,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_x,tol_p,tol_of,second_order)
    elseif  !isnothing(p0) && !isnothing(x0)
        p0,_ = promote(p0,first(x0))
        T = eltype(p0)
        x0 = convert(Vector{T},x0)
        return FugDewPressure{T}(vol0,p0,x0,noncondensables,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_x,tol_p,tol_of,second_order)
    else
        throw(error("invalid specification for dew pressure"))
    end
end

function dew_pressure_impl(model::RestrictedEquilibriaModel,T,y,method::FugDewPressure)
    wrapper = PTFlashWrapper(model,NaN,T,y,:vle)
    return dew_pressure_impl(wrapper,T,y,method)
end

function dew_pressure_impl(model::CompositeModel,T,y,method::FugDewPressure)
    return dew_pressure_impl(model.fluid,T,y,method)
end

function dew_pressure_impl(model::EoSModel, T, y ,method::FugDewPressure)
    noncondensables = method.noncondensables
    condensables = comps_in_equilibria(component_list(model),noncondensables)
    _vol0,_p0,_x0 = method.vol0,method.p0,method.x0
    p0,volx,voly,x0 = dew_pressure_init(model,T,y,_vol0,_p0,_x0,condensables)
    vol0 = (volx,voly)
    model_x,_ = index_reduction(model,condensables)
    x0 = x0[condensables]
    x0 = x0/sum(x0)

    data = FugData(FugEnum.DEW_PRESSURE,
                method.itmax_ss,
                method.itmax_newton,
                method.tol_p,
                method.tol_x,
                method.tol_of,
                method.second_order,
                false)
    
    cache = fug_bubbledew_cache(model_x,model,p0,T,x0,y,data)

    if all(condensables)
        converged,res = _fug_OF_ss(model,p0,T,x0,y,vol0,data,cache)
    else
        converged,res = _fug_OF_ss(model_x,model,p0,T,x0,y,vol0,condensables,data,cache)
    end

    p,T,x,y,vol = res
    volx,voly = vol

    if converged
        return p,volx,voly,index_expansion(x,condensables)
    end

    lnK,K,w,w_old,w_calc,w_restart,vol_cache,Hϕx,Hϕy = cache
    inc0 = vcat(lnK, log(p))
    vol_cache[] = (volx,voly)
    opts = NLSolvers.NEqOptions(method)
    if all(condensables)
        problem = _fug_OF_neq(model,T,y,data,cache)
        sol = Solvers.nlsolve(problem, inc0, Solvers.LineSearch(Solvers.Newton2(inc0)),opts)
        inc = Solvers.x_sol(sol)
        !all(<(sol.options.f_abstol),sol.info.best_residual) && (inc .= NaN)
    else
        problem = _fug_OF_neq(model_x,model,T,y,condensables,data,cache)
        sol = Solvers.nlsolve(problem, inc0, Solvers.LineSearch(Solvers.Newton2(inc0)),opts)
        inc = Solvers.x_sol(sol)
        !all(<(sol.options.f_abstol),sol.info.best_residual) && (inc .= NaN)
    end

    lnp = inc[end]
    lnK = @view inc[1:(end-1)]
    x_r = @view(y[condensables]) ./ exp.(lnK)
    x = index_expansion(x_r,condensables)
    p = exp(lnp)
    volx,voly = vol_cache[]

     return p, volx, voly, x
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
    second_order::Bool
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
    tol_of = 1e-8,
    second_order = true)

    if T0 == x0 == vol0 == nothing
        return FugDewTemperature{Nothing}(vol0,T0,x0,noncondensables,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_x,tol_T,tol_of,second_order)
    elseif (T0 == x0 == nothing) && !isnothing(vol0)
        vl,vv = promote(vol0[1],vol0[2])
        return FugDewTemperature{typeof(vl)}(vol0,T0,x0,noncondensables,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_x,tol_T,tol_of,second_order)
    elseif (vol0 == x0 == nothing) && !isnothing(T0)
        T0 = float(T0)
        return FugDewTemperature{typeof(T0)}(vol0,T0,x0,noncondensables,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_x,tol_T,tol_of,second_order)
    elseif (T0 == vol0 == nothing) && !isnothing(x0)
        T = eltype(x0)
        return FugDewTemperature{T}(vol0,T0,x0,noncondensables,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_x,tol_T,tol_of,second_order)
    elseif !isnothing(vol0) && !isnothing(T0) && !isnothing(x0)
        vl,vv,T0,_ = promote(vol0[1],vol0[2],T0,first(x0))
        T = eltype(vl)
        x0 = convert(Vector{T},x0)
        return FugDewTemperature{T}(vol0,T0,x0,noncondensables,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_x,tol_T,tol_of,second_order)
    elseif !isnothing(vol0) && !isnothing(x0)
        vl,vv,_ = promote(vol0[1],vol0[2],first(x0))
        T = eltype(vl)
        x0 = convert(Vector{T},x0)
        return FugDewTemperature{T}(vol0,T0,x0,noncondensables,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_x,tol_T,tol_of,second_order)
    elseif  !isnothing(T0) && !isnothing(x0)
        T0,_ = promote(T0,first(x0))
        T = eltype(T0)
        x0 = convert(Vector{T},x0)
        return FugDewTemperature{T}(vol0,T0,x0,noncondensables,f_limit,atol,rtol,max_iters,itmax_newton,itmax_ss,tol_x,tol_T,tol_of,second_order)
    else
        throw(error("invalid specification for bubble temperature"))
    end
end

function dew_temperature_impl(model::RestrictedEquilibriaModel,p,y,method::FugDewTemperature)
    wrapper = PTFlashWrapper(model,p,NaN,y,:vle)
    return dew_temperature_impl(wrapper,p,y,method)
end

function dew_temperature_impl(model::CompositeModel,p,y,method::FugDewTemperature)
    return dew_temperature_impl(model.fluid,p,y,method)
end

function dew_temperature_impl(model::EoSModel, p, y, method::FugDewTemperature)
    noncondensables = method.noncondensables
    condensables = comps_in_equilibria(component_list(model),noncondensables)
    _vol0,_T0,_x0 = method.vol0,method.T0,method.x0
    T0,volx,voly,x0 = dew_temperature_init(model,p,y,_vol0,_T0,_x0,condensables)
    model_x,_ = index_reduction(model,condensables)
    vol0 = (volx,voly)
    x0 = x0[condensables]
    x0 ./= sum(x0)

    data = FugData(FugEnum.DEW_TEMPERATURE,
                method.itmax_ss,
                method.itmax_newton,
                method.tol_T,
                method.tol_x,
                method.tol_of,
                method.second_order,
                false)
    
    cache = fug_bubbledew_cache(model_x,model,p,T0,x0,y,data)

    if all(condensables)
        converged,res = _fug_OF_ss(model,p,T0,x0,y,vol0,data,cache)
    else
        converged,res = _fug_OF_ss(model_x,model,p,T0,x0,y,vol0,condensables,data,cache)
    end
    
    p,T,x,y,vol = res
    volx,voly = vol

    if converged
        return T,volx,voly,index_expansion(x,condensables)
    end

    lnK,K,w,w_old,w_calc,w_restart,vol_cache,Hϕx,Hϕy = cache
    inc0 = vcat(lnK, log(T))
    vol_cache[] = (volx,voly)
    opts = NLSolvers.NEqOptions(method)
    if all(condensables)
        problem = _fug_OF_neq(model,p,y,data,cache)
        sol = Solvers.nlsolve(problem, inc0, Solvers.LineSearch(Solvers.Newton2(inc0)),opts)
        inc = Solvers.x_sol(sol)
        !all(<(sol.options.f_abstol),sol.info.best_residual) && (inc .= NaN)
    else
        problem = _fug_OF_neq(model_x,model,p,y,condensables,data,cache)
        sol = Solvers.nlsolve(problem, inc0, Solvers.LineSearch(Solvers.Newton2(inc0)),opts)
        inc = Solvers.x_sol(sol)
        !all(<(sol.options.f_abstol),sol.info.best_residual) && (inc .= NaN)
    end

    lnT = inc[end]
    lnK = @view inc[1:(end-1)]
    x_r = @view(y[condensables]) ./ exp.(lnK)
    x = index_expansion(x_r,condensables)
    T = exp(lnT)
    volx,voly = vol_cache[]
    
    return T, volx, voly, x
end

export FugDewPressure, FugDewTemperature