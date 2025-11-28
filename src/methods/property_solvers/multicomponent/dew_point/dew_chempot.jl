"""
    ChemPotDewPressure(kwargs...)

Function to compute [`dew_pressure`](@ref) via chemical potentials.
It directly solves the equality of chemical potentials system of equations.

Inputs:
- `x0 = nothing`: optional, initial guess for the liquid phase composition
- `p0 = nothing`: optional, initial guess for the dew pressure `[Pa]`
- `vol0 = nothing`: optional, initial guesses for the liquid and vapor phase volumes `[m³]`
- `atol = 1e-8`: optional, absolute tolerance of the non linear system of equations
- `rtol = 1e-12`: optional, relative tolerance of the non linear system of equations
- `max_iters = 1000`: optional, maximum number of iterations
- `noncondensables = nothing`: optional, Vector of strings containing non condensable compounds. those will be set to zero on the liquid phase.
"""
struct ChemPotDewPressure{T} <: DewPointMethod
    vol0::Union{Nothing,Tuple{T,T}}
    p0::Union{Nothing,T}
    x0::Union{Nothing,Vector{T}}
    noncondensables::Union{Nothing,Vector{String}}
    f_limit::Float64
    atol::Float64
    rtol::Float64
    max_iters::Int
    ss::Bool
end

function ChemPotDewPressure(;vol0 = nothing,
                                p0 = nothing,
                                x0 = nothing,
                                noncondensables = nothing,
                                f_limit = 0.0,
                                atol = 1e-8,
                                rtol = 1e-12,
                                max_iters = 10^3,
                                ss = false)

    if p0 == x0 == vol0 == nothing
        return ChemPotDewPressure{Nothing}(vol0,p0,x0,noncondensables,f_limit,atol,rtol,max_iters,ss)
    elseif (p0 == x0 == nothing) && !isnothing(vol0)
        vl,vv = promote(vol0[1],vol0[2])
        return ChemPotDewPressure{typeof(vl)}(vol0,p0,x0,noncondensables,f_limit,atol,rtol,max_iters,ss)
    elseif (vol0 == x0 == nothing) && !isnothing(p0)
        p0 = float(p0)
        return ChemPotDewPressure{typeof(p0)}(vol0,p0,x0,noncondensables,f_limit,atol,rtol,max_iters,ss)
    elseif (p0 == vol0 == nothing) && !isnothing(x0)
        T = eltype(x0)
        return ChemPotDewPressure{T}(vol0,p0,x0,noncondensables,f_limit,atol,rtol,max_iters,ss)
    elseif !isnothing(vol0) && !isnothing(p0) && !isnothing(x0)
        vl,vv,p0,_ = promote(vol0[1],vol0[2],p0,first(x0))
        T = eltype(vl)
        x0 = convert(Vector{T},x0)
        return ChemPotDewPressure{T}(vol0,p0,x0,noncondensables,f_limit,atol,rtol,max_iters,ss)
    elseif !isnothing(vol0) && !isnothing(x0)
        vl,vv,_ = promote(vol0[1],vol0[2],first(x0))
        T = eltype(vl)
        x0 = convert(Vector{T},x0)
        return ChemPotDewPressure{T}(vol0,p0,x0,noncondensables,f_limit,atol,rtol,max_iters,ss)
    elseif  !isnothing(p0) && !isnothing(x0)
        p0,_ = promote(p0,first(x0))
        T = eltype(p0)
        x0 = convert(Vector{T},x0)
        return ChemPotDewPressure{T}(vol0,p0,x0,noncondensables,f_limit,atol,rtol,max_iters,ss)
    else
        throw(error("invalid specification for dew pressure"))
    end
end

function dew_pressure_impl(model::EoSModel, T, y, method::ChemPotDewPressure)
    is_non_condensable = !isnothing(method.noncondensables)
    condensables = comps_in_equilibria(component_list(model),method.noncondensables)
    model_x,_ = index_reduction(model,condensables)
    p0,vl0,vv0,x00 = dew_pressure_init(model,T,y,method.vol0,method.p0,method.x0,condensables)
    x0 = x00[condensables]
    data = FugEnum.DEW_PRESSURE
    
    neq = count(condensables)
    cache = Clapeyron.fug_bubbledew_cache(model_x,model,T,T,y,y,Val{false}())
    x_r,_,_,_,_,_,p_cache,_,_ = cache
    inc0 = similar(x_r,neq+2)
    inc0[1:neq] .= log.(@view(y[condensables]) ./ x0)
    inc0[neq+1] = log(vl0)
    inc0[neq+2] = log(vv0)
    opts = NLSolvers.NEqOptions(method)
    if !is_non_condensable
        problem = _mu_OF_neq(model,T,y,data,cache)
        sol = Solvers.nlsolve(problem, inc0, Solvers.LineSearch(Solvers.Newton2(inc0),Static(1.0)),opts)
        inc = Solvers.x_sol(sol)
        !all(<(sol.options.f_abstol),sol.info.best_residual) && (inc .= NaN)
    else
        problem = _mu_OF_neq(model_x,model,T,y,condensables,data,cache)
        sol_vol = Solvers.nlsolve(problem, inc0, Solvers.LineSearch(Solvers.Newton2(inc0),Static(1.0)),opts)
        inc = Solvers.x_sol(sol_vol)
        !all(<(sol_vol.options.f_abstol),sol_vol.info.best_residual) && (inc .= NaN)
    end
    lnK = @view inc[1:neq]
    x_r .= @view(y[condensables]) ./  exp.(lnK)
    x = index_expansion(x_r,condensables)
    vl = exp(inc[neq+1])
    vv = exp(inc[neq+2])
    px,py = p_cache[]
    p = 0.5*(px + py)
    return (p, vl, vv, x)
end

function Obj_dew_pressure(model::EoSModel,model_x, F, T, ηl, ηv, x, y, _view, xx_i)
    xx = FractionVector(x,xx_i)
    vl = v_from_η(model_x,ηl,T,xx)
    vv = v_from_η(model,ηv,T,y)
    v = (vv,vl)
    w = (y,xx)
    if all(_view)
        return μp_equality2(model, nothing, F, Tspec(T), v, w, _view)
    else
        return μp_equality2(model, model_x, F, Tspec(T), v, w, _view)
    end
end

"""
    ChemPotDewTemperature(kwargs...)

Function to compute [`dew_temperature`](@ref) via chemical potentials.
It directly solves the equality of chemical potentials system of equations.

Inputs:
- `x0 = nothing`: optional, initial guess for the liquid phase composition
- `T0  =nothing`: optional, initial guess for the dew temperature `[K]`
- `vol0 = nothing`: optional, initial guesses for the liquid and vapor phase volumes `[m³]`
- `atol = 1e-8`: optional, absolute tolerance of the non linear system of equations
- `rtol = 1e-12`: optional, relative tolerance of the non linear system of equations
- `max_iters = 1000`: optional, maximum number of iterations
- `noncondensables = nothing`: optional, Vector of strings containing non condensable compounds. those will be set to zero on the liquid phase.
"""
struct ChemPotDewTemperature{T} <: DewPointMethod
    vol0::Union{Nothing,Tuple{T,T}}
    T0::Union{Nothing,T}
    x0::Union{Nothing,Vector{T}}
    noncondensables::Union{Nothing,Vector{String}}
    f_limit::Float64
    atol::Float64
    rtol::Float64
    max_iters::Int
    ss::Bool
end

function ChemPotDewTemperature(;vol0 = nothing,
    T0 = nothing,
    x0 = nothing,
    noncondensables = nothing,
    f_limit = 0.0,
    atol = 1e-8,
    rtol = 1e-12,
    max_iters = 10^3,
    ss = false)

    if T0 == x0 == vol0 == nothing
        return ChemPotDewTemperature{Nothing}(vol0,T0,x0,noncondensables,f_limit,atol,rtol,max_iters,ss)
    elseif (T0 == x0 == nothing) && !isnothing(vol0)
        vl,vv = promote(vol0[1],vol0[2])
        return ChemPotDewTemperature{typeof(vl)}(vol0,T0,x0,noncondensables,f_limit,atol,rtol,max_iters,ss)
    elseif (vol0 == x0 == nothing) && !isnothing(T0)
        T0 = float(T0)
        return ChemPotDewTemperature{typeof(T0)}(vol0,T0,x0,noncondensables,f_limit,atol,rtol,max_iters,ss)
    elseif (T0 == vol0 == nothing) && !isnothing(x0)
        T = eltype(x0)
        return ChemPotDewTemperature{T}(vol0,T0,x0,noncondensables,f_limit,atol,rtol,max_iters,ss)
    elseif !isnothing(vol0) && !isnothing(T0) && !isnothing(x0)
        vl,vv,T0,_ = promote(vol0[1],vol0[2],T0,first(x0))
        T = eltype(vl)
        x0 = convert(Vector{T},x0)
        return ChemPotDewTemperature{T}(vol0,T0,x0,noncondensables,f_limit,atol,rtol,max_iters,ss)
    elseif !isnothing(vol0) && !isnothing(x0)
        vl,vv,_ = promote(vol0[1],vol0[2],first(x0))
        T = eltype(vl)
        x0 = convert(Vector{T},x0)
        return ChemPotDewTemperature{T}(vol0,T0,x0,noncondensables,f_limit,atol,rtol,max_iters,ss)
    elseif  !isnothing(T0) && !isnothing(x0)
        T0,_ = promote(T0,first(x0))
        T = eltype(T0)
        x0 = convert(Vector{T},x0)
        return ChemPotDewTemperature{T}(vol0,T0,x0,noncondensables,f_limit,atol,rtol,max_iters,ss)
    else
        throw(error("invalid specification for bubble temperature"))
    end
end

function dew_temperature_impl(model::EoSModel, p, y, method::ChemPotDewTemperature)
    is_non_condensable = !isnothing(method.noncondensables)
    condensables = comps_in_equilibria(component_list(model),method.noncondensables)
    model_x,_ = index_reduction(model,condensables)
    T0,vl0,vv0,x00 = dew_temperature_init(model,p,y,method.vol0,method.T0,method.x0,condensables)
    x0 = x00[condensables]
    data = FugEnum.DEW_TEMPERATURE
    
    neq = count(condensables)
    cache = Clapeyron.fug_bubbledew_cache(model_x,model,p,p,y,y,Val{true}())
    x_r,_,_,_,_,_,_,_,_ = cache
    inc0 = similar(x_r,neq+3)
    inc0[1:neq] .= log.(@view(y[condensables]) ./ x0)
    inc0[neq+1] = log(vl0)
    inc0[neq+2] = log(vv0)
    inc0[neq+3] = log(T0)
    opts = NLSolvers.NEqOptions(method)
    if !is_non_condensable
        problem = _mu_OF_neq(model,p,y,data,cache)
        sol = Solvers.nlsolve(problem, inc0, Solvers.LineSearch(Solvers.Newton2(inc0),Static(1.0)),opts)
        inc = Solvers.x_sol(sol)
        !all(<(sol.options.f_abstol),sol.info.best_residual) && (inc .= NaN)
    else
        problem = _mu_OF_neq(model_x,model,p,y,condensables,data,cache)
        sol_vol = Solvers.nlsolve(problem, inc0, Solvers.LineSearch(Solvers.Newton2(inc0),Static(1.0)),opts)
        inc = Solvers.x_sol(sol_vol)
        !all(<(sol_vol.options.f_abstol),sol_vol.info.best_residual) && (inc .= NaN)
    end
    lnK = @view inc[1:neq]
    x_r .= @view(y[condensables]) ./  exp.(lnK)
    x = index_expansion(x_r,condensables)
    vl = exp(inc[neq+1])
    vv = exp(inc[neq+2])
    T = exp(inc[neq+3])
    return (T, vl, vv, x)
end

function Obj_dew_temperature(model::EoSModel,model_x, F, p, T, ηl, ηv, x, y, _view,xx_i)
    vv = v_from_η(model,ηv,T,y)
    xx = FractionVector(x,xx_i)
    w = (y,xx)
    vl = v_from_η(model_x,ηl,T,xx)
    v = (vv,vl)
    if all(_view)
        return μp_equality2(model, nothing, F, Pspec(p,T), v, w, _view)
    else
        return μp_equality2(model, model_x, F, Pspec(p,T), v, w, _view)
    end
end

export ChemPotDewPressure, ChemPotDewTemperature
