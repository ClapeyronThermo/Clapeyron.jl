"""
    ChemPotDewPressure(kwargs...)

Function to compute [`dew_pressure`](@ref) via chemical potentials.
It directly solves the equality of chemical potentials system of equations.

Inputs:
- `x0 = nothing`: optional, initial guess for the liquid phase composition
- `p0 = nothing`: optional, initial guess for the dew pressure [`Pa`]
- `vol0 = nothing`: optional, initial guesses for the liquid and vapor phase volumes
- `atol = 1e-8`: optional, absolute tolerance of the non linear system of equations
- `rtol = 1e-12`: optional, relative tolerance of the non linear system of equations
- `max_iters = 10000`: optional, maximum number of iterations
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
end

function ChemPotDewPressure(;vol0 = nothing,
                                p0 = nothing,
                                x0 = nothing,
                                noncondensables = nothing,
                                f_limit = 0.0,
                                atol = 1e-8,
                                rtol = 1e-12,
                                max_iters = 10^4)

    if p0 == x0 == vol0 == nothing
        return ChemPotDewPressure{Nothing}(vol0,p0,x0,noncondensables,f_limit,atol,rtol,max_iters)
    elseif (p0 == x0 == nothing) && !isnothing(vol0)
        vl,vv = promote(vol0[1],vol0[2])
        return ChemPotDewPressure{typeof(vl)}(vol0,p0,x0,noncondensables,f_limit,atol,rtol,max_iters)
    elseif (vol0 == x0 == nothing) && !isnothing(p0)
        p0 = float(p0)
        return ChemPotDewPressure{typeof(p0)}(vol0,p0,x0,noncondensables,f_limit,atol,rtol,max_iters)
    elseif (p0 == vol0 == nothing) && !isnothing(x0)
        T = eltype(x0)
        return ChemPotDewPressure{T}(vol0,p0,x0,noncondensables,f_limit,atol,rtol,max_iters)
    elseif !isnothing(vol0) && !isnothing(p0) && !isnothing(x0)
        vl,vv,p0,_ = promote(vol0[1],vol0[2],p0,first(x0))
        T = eltype(vl)
        x0 = convert(Vector{T},x0)
        return ChemPotDewPressure{T}(vol0,p0,x0,noncondensables,f_limit,atol,rtol,max_iters)
    elseif !isnothing(vol0) && !isnothing(x0)
        vl,vv,_ = promote(vol0[1],vol0[2],first(x0))
        T = eltype(vl)
        x0 = convert(Vector{T},x0)
        return ChemPotDewPressure{T}(vol0,p0,x0,noncondensables,f_limit,atol,rtol,max_iters)
    elseif  !isnothing(p0) && !isnothing(x0)
        p0,_ = promote(p0,first(x0))
        T = eltype(p0)
        x0 = convert(Vector{T},x0)
        return ChemPotDewPressure{T}(vol0,p0,x0,noncondensables,f_limit,atol,rtol,max_iters)
    else
        throw(error("invalid specification for dew pressure"))
    end
end

function dew_pressure_impl(model::EoSModel, T, y,method::ChemPotDewPressure)

    if !isnothing(method.noncondensables)
        condensables = [!in(x,method.noncondensables) for x in model.components]
    else
        condensables = fill(true,length(model))
    end

    _vol0,_p0,_x0 = method.vol0,method.p0,method.x0
    p0,vl,vv,x0 = dew_pressure_init(model,T,y,_vol0,_p0,_x0,condensables)

    if !isnothing(method.noncondensables)
        model_x,condensables = index_reduction(model,condensables)
        x0 = x0[condensables]
    else
        model_x = nothing
    end
    v0 = vcat(log10(vl),log10(vv),x0[1:end-1])
    pmix = p_scale(model,y)
    f!(F,z) = Obj_dew_pressure(model,model_x, F, T, exp10(z[1]), exp10(z[2]), z[3:end],y,pmix,condensables)
    r  =Solvers.nlsolve(f!,v0,LineSearch(Newton()),NLSolvers.NEqOptions(method))
    sol = Solvers.x_sol(r)
    v_l = exp10(sol[1])
    v_v = exp10(sol[2])
    x_r = FractionVector(sol[3:end])
    P_sat = pressure(model,v_v,T,y)
    x = index_expansion(collect(x_r),condensables)
    return (P_sat, v_l, v_v, x)
end

function Obj_dew_pressure(model::EoSModel,model_x, F, T, v_l, v_v, x, y,ps,_view)
    return μp_equality(model,model_x, F, T, v_v, v_l, y,FractionVector(x),ps,_view)
end

"""
    ChemPotDewTemperature(kwargs...)

Function to compute [`dew_temperature`](@ref) via chemical potentials.
It directly solves the equality of chemical potentials system of equations.

Inputs:
- `x0 = nothing`: optional, initial guess for the liquid phase composition
- `T0  =nothing`: optional, initial guess for the dew temperature [`K`]
- `vol0 = nothing`: optional, initial guesses for the liquid and vapor phase volumes
- `atol = 1e-8`: optional, absolute tolerance of the non linear system of equations
- `rtol = 1e-12`: optional, relative tolerance of the non linear system of equations
- `max_iters = 10000`: optional, maximum number of iterations
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
end

function ChemPotDewTemperature(;vol0 = nothing,
    T0 = nothing,
    x0 = nothing,
    noncondensables = nothing,
    f_limit = 0.0,
    atol = 1e-8,
    rtol = 1e-12,
    max_iters = 10^4)

    if T0 == x0 == vol0 == nothing
        return ChemPotDewTemperature{Nothing}(vol0,T0,x0,noncondensables,f_limit,atol,rtol,max_iters)
    elseif (T0 == x0 == nothing) && !isnothing(vol0)
        vl,vv = promote(vol0[1],vol0[2])
        return ChemPotDewTemperature{typeof(vl)}(vol0,T0,x0,noncondensables,f_limit,atol,rtol,max_iters)
    elseif (vol0 == x0 == nothing) && !isnothing(T0)
        T0 = float(T0)
        return ChemPotDewTemperature{typeof(T0)}(vol0,T0,x0,noncondensables,f_limit,atol,rtol,max_iters)
    elseif (T0 == vol0 == nothing) && !isnothing(x0)
        T = eltype(x0)
        return ChemPotDewTemperature{T}(vol0,T0,x0,noncondensables,f_limit,atol,rtol,max_iters)
    elseif !isnothing(vol0) && !isnothing(T0) && !isnothing(x0)
        vl,vv,T0,_ = promote(vol0[1],vol0[2],T0,first(x0))
        T = eltype(vl)
        x0 = convert(Vector{T},x0)
        return ChemPotDewTemperature{T}(vol0,T0,x0,noncondensables,f_limit,atol,rtol,max_iters)
    elseif !isnothing(vol0) && !isnothing(x0)
        vl,vv,_ = promote(vol0[1],vol0[2],first(x0))
        T = eltype(vl)
        x0 = convert(Vector{T},x0)
        return ChemPotDewTemperature{T}(vol0,T0,x0,noncondensables,f_limit,atol,rtol,max_iters)
    elseif  !isnothing(T0) && !isnothing(x0)
        T0,_ = promote(T0,first(x0))
        T = eltype(T0)
        x0 = convert(Vector{T},x0)
        return ChemPotDewTemperature{T}(vol0,T0,x0,noncondensables,f_limit,atol,rtol,max_iters)
    else
        throw(error("invalid specification for bubble temperature"))
    end
end

function dew_temperature_impl(model::EoSModel,p,y,method::ChemPotDewTemperature)
    T0,vl,vv,x0 = dew_temperature_init(model,p,y,method.vol0,method.T0,method.x0)
    if !isnothing(method.noncondensables)
        condensables = [!in(x,method.noncondensables) for x in model.components]
        model_x,condensables = index_reduction(model,condensables)
        x0 = x0[condensables]
    else
        condensables = fill(true,length(model))
        model_x = nothing
    end
    v0 = vcat(T0,log10(vl),log10(vv),x0[1:end-1])
    pmix = p_scale(model,y)
    f!(F,z) = Obj_dew_temperature(model,model_x, F, p, z[1], exp10(z[2]), exp10(z[3]), z[4:end],y, pmix, condensables)
    r  =Solvers.nlsolve(f!,v0,LineSearch(Newton()),NLSolvers.NEqOptions(method))
    sol = Solvers.x_sol(r)
    T   = sol[1]
    v_l = exp10(sol[2])
    v_v = exp10(sol[3])
    x_r = FractionVector(sol[4:end])
    x = index_expansion(x_r,condensables)
    return T, v_l, v_v, x
end

function Obj_dew_temperature(model::EoSModel,model_x, F, p, T, v_l, v_v, x, y, ps, _view)
    Ts = T_scale(model,y)
    F = μp_equality(model, model_x, F, T, v_v, v_l, y, FractionVector(x), ps, _view, Ts)
    F[end] = (pressure(model,v_v,T,y) - p)/ps
    return F
end

export ChemPotDewPressure, ChemPotDewTemperature