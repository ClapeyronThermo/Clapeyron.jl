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

function dew_pressure_impl(model::EoSModel, T, y,method::ChemPotDewPressure)

    is_non_condensable = !isnothing(method.noncondensables)
    condensables = comps_in_equilibria(model.components,method.noncondensables)
    model_x,_ = index_reduction(model,condensables)
    p0,vl,vv,x0 = dew_pressure_init(model,T,y,method.vol0,method.p0,method.x0,condensables)

    if is_non_condensable
        ηl0 = η_from_v(model_x,vl,T,x0)
    else
        ηl0 = η_from_v(model,vl,T,x0)
    end

    ηv0 = η_from_v(model,vv,T,y)
    _,idx_max = findmax(x0)
    v0 = vcat(ηl0,ηv0,deleteat(x0,idx_max)) #select component with highest fraction as pivot
    f!(F,z) = Obj_dew_pressure(model,model_x, F, T, z[1], z[2], z[3:end], y, condensables, idx_max)
    r = Solvers.nlsolve(f!,v0,LineSearch(Newton2(v0)),NLSolvers.NEqOptions(method))
    sol = Solvers.x_sol(r)
    !all(<(r.options.f_abstol),r.info.best_residual) && (sol .= NaN)
    x_r = FractionVector(sol[3:end],idx_max)
    v_l = v_from_η(model,model_x,sol[1],T,x_r)
    v_v = v_from_η(model,sol[2],T,y)
    P_sat = pressure(model,v_v,T,y)
    x = index_expansion(x_r,condensables)
    return (P_sat, v_l, v_v, x)
end

function Obj_dew_pressure(model::EoSModel,model_x, F, T, ηl, ηv, x, y, _view, xx_i)
    xx = FractionVector(x,xx_i)
    vl = v_from_η(model,model_x,ηl,T,xx)
    vv = v_from_η(model,ηv,T,y)
    v = (vv,vl)
    w = (y,xx)
    if all(_view)
        return μp_equality2(model,nothing, F, Tspec(T), v, w, _view)
    else
        return μp_equality2(model,model_x, F, Tspec(T), v, w, _view)
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

function dew_temperature_impl(model::EoSModel,p,y,method::ChemPotDewTemperature)

    is_non_condensable = !isnothing(method.noncondensables)
    condensables = comps_in_equilibria(model.components,method.noncondensables)
    model_x,_ = index_reduction(model,condensables)
    T0,vl,vv,x0 = dew_temperature_init(model,T,y,method.vol0,method.T0,method.x0,condensables)

    if is_non_condensable
        ηl0 = η_from_v(model_x,vl,T0,x0)
    else
        ηl0 = η_from_v(model,vl,T0,x0)
    end

    ηl = η_from_v(model,model_x,vl,T0,x0)
    ηv = η_from_v(model,vv,T0,y)
    _,idx_max = findmax(x0)
    v0 = vcat(T0,ηl,ηv,deleteat(x0,idx_max)) #select component with highest fraction as pivot
    f!(F,z) = Obj_dew_temperature(model,model_x, F, p, z[1], z[2], z[3], z[4:end], y, condensables,idx_max)
    r = Solvers.nlsolve(f!,v0,LineSearch(Newton2(v0)),NLSolvers.NEqOptions(method))
    sol = Solvers.x_sol(r)
    !all(<(r.options.f_abstol),r.info.best_residual) && (sol .= NaN)
    T   = sol[1]
    x_r = FractionVector(sol[4:end],idx_max)
    v_l = v_from_η(model,model_x,sol[2],T,x_r)
    v_v = v_from_η(model,sol[3],T,y)
    x = index_expansion(x_r,condensables)
    return T, v_l, v_v, x
end

function Obj_dew_temperature(model::EoSModel,model_x, F, p, T, ηl, ηv, x, y, _view,xx_i)
    vv = v_from_η(model,ηv,T,y)
    xx = FractionVector(x,xx_i)
    w = (y,xx)
    vl = v_from_η(model,model_x,ηl,T,xx)
    v = (vv,vl)
    if all(_view)
        return μp_equality2(model, nothing, F, Pspec(p,T), v, w, _view)
    else
        return μp_equality2(model, model_x, F, Pspec(p,T), v, w, _view)
    end
end

export ChemPotDewPressure, ChemPotDewTemperature
