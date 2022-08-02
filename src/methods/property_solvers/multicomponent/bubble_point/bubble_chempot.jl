## Bubble pressure solver
"""
    ChemPotBubblePressure(kwargs...)  

Function to compute [`bubble_pressure`](@ref) via chemical potentials. 
It directly solves the equality of chemical potentials system of equations.

Inputs:
- `y0 = nothing`: optional, initial guess for the vapor phase composition
- `p0 = nothing`: optional, initial guess for the bubble pressure [`Pa`]
- `vol0 = nothing`: optional, initial guesses for the liquid and vapor phase volumes
- `atol = 1e-8`: optional, absolute tolerance of the non linear system of equations
- `rtol = 1e-12`: optional, relative tolerance of the non linear system of equations
- `max_iters = 10000`: optional, maximum number of iterations
- `nonvolatiles = nothing`: optional, Vector of strings containing non volatile compounds. those will be set to zero on the vapour phase.
"""
struct ChemPotBubblePressure{T} <: BubblePointMethod
    vol0::Union{Nothing,Tuple{T,T}}
    p0::Union{Nothing,T}
    y0::Union{Nothing,Vector{T}}
    nonvolatiles::Union{Nothing,Vector{String}}
    f_limit::Float64
    atol::Float64
    rtol::Float64
    max_iters::Int
end

function ChemPotBubblePressure(;vol0 = nothing,
                                p0 = nothing,
                                y0 = nothing,
                                nonvolatiles = nothing,
                                f_limit = 0.0,
                                atol = 1e-8,
                                rtol = 1e-12,
                                max_iters = 10^4)

    if p0 == y0 == vol0 == nothing
        return ChemPotBubblePressure{Nothing}(vol0,p0,y0,nonvolatiles,f_limit,atol,rtol,max_iters)
    elseif (p0 == y0 == nothing) && !isnothing(vol0)
        vl,vv = promote(vol0[1],vol0[2])
        return ChemPotBubblePressure{typeof(vl)}(vol0,p0,y0,nonvolatiles,f_limit,atol,rtol,max_iters)
    elseif (vol0 == y0 == nothing) && !isnothing(p0)
        p0 = float(p0)
        return ChemPotBubblePressure{typeof(p0)}(vol0,p0,y0,nonvolatiles,f_limit,atol,rtol,max_iters)
    elseif (p0 == vol0 == nothing) && !isnothing(y0)
        T = eltype(y0)
        return ChemPotBubblePressure{T}(vol0,p0,y0,nonvolatiles,f_limit,atol,rtol,max_iters)
    elseif !isnothing(vol0) && !isnothing(p0) && !isnothing(y0)
        vl,vv,p0,_ = promote(vol0[1],vol0[2],p0,first(y0))
        T = eltype(vl)
        y0 = convert(Vector{T},y0)
        return ChemPotBubblePressure{T}(vol0,p0,y0,nonvolatiles,f_limit,atol,rtol,max_iters)
    elseif !isnothing(vol0) && !isnothing(y0)
        vl,vv,_ = promote(vol0[1],vol0[2],first(y0))
        T = eltype(vl)
        y0 = convert(Vector{T},y0)
        return ChemPotBubblePressure{T}(vol0,p0,y0,nonvolatiles,f_limit,atol,rtol,max_iters)
    elseif  !isnothing(p0) && !isnothing(y0)
        p0,_ = promote(p0,first(y0))
        T = eltype(p0)
        y0 = convert(Vector{T},y0)
        return ChemPotBubblePressure{T}(vol0,p0,y0,nonvolatiles,f_limit,atol,rtol,max_iters)
    else
        throw(error("invalid specification for bubble pressure"))
    end
end

function bubble_pressure_impl(model::EoSModel, T, x,method::ChemPotBubblePressure)
    p0,vl,vv,y0 = bubble_pressure_init(model,T,x,method.vol0,method.p0,method.y0)
    if !isnothing(method.nonvolatiles)
        volatiles = [!in(x,method.nonvolatiles) for x in model.components]
        model_y,volatiles = index_reduction(model,volatiles)
        y0 = y0[volatiles]
        ts = T_scales(model_y,y0)
    else
        volatiles = fill(true,length(model))
        model_y = nothing
        ts = T_scales(model)
    end
    v0 = vcat(log10(vl),log10(vv),y0[1:end-1])    
    pmix = p_scale(model,x)
    f!(F,z) = Obj_bubble_pressure(model,model_y, F, T, exp10(z[1]),exp10(z[2]),x,z[3:end],ts,pmix,volatiles)
    r  =Solvers.nlsolve(f!,v0,LineSearch(Newton()),NLSolvers.NEqOptions(method))
    sol = Solvers.x_sol(r)
    v_l = exp10(sol[1])
    v_v = exp10(sol[2])
    y = FractionVector(sol[3:end])
    y = index_expansion(collect(y),volatiles)
    P_sat = pressure(model,v_l,T,x)
    return (P_sat, v_l, v_v, y)
end


function Obj_bubble_pressure(model::EoSModel, model_y, F, T, v_l, v_v, x, y,ts,ps,_view)
    return μp_equality(model,model_y ,F, T, v_l, v_v, x, FractionVector(y),ts,ps,_view)
end

#used by LLE_pressure
function Obj_bubble_pressure(model::EoSModel, F, T, v_l, v_v, x, y,ts,ps)
    return Obj_bubble_pressure(model, nothing, F, T, v_l, v_v, x, y,ts,ps,nothing)
end


struct ChemPotBubbleTemperature{T} <: BubblePointMethod
    vol0::Union{Nothing,Tuple{T,T}}
    T0::Union{Nothing,T}
    y0::Union{Nothing,Vector{T}}
    nonvolatiles::Union{Nothing,Vector{String}}
    f_limit::Float64
    atol::Float64
    rtol::Float64
    max_iters::Int
end

"""
    ChemPotBubbleTemperature(kwargs...)  

Function to compute [`bubble_temperature`](@ref) via chemical potentials. 
It directly solves the equality of chemical potentials system of equations.

Inputs:
- `y = nothing`: optional, initial guess for the vapor phase composition.
- `T0 = nothing`: optional, initial guess for the bubble temperature [`K`].
- `vol0 = nothing`: optional, initial guesses for the liquid and vapor phase volumes
- `atol = 1e-8`: optional, absolute tolerance of the non linear system of equations
- `rtol = 1e-12`: optional, relative tolerance of the non linear system of equations
- `max_iters = 10000`: optional, maximum number of iterations
- `nonvolatiles = nothing`: optional, Vector of strings containing non volatile compounds. those will be set to zero on the vapour phase.
"""
function ChemPotBubbleTemperature(;vol0 = nothing,
    T0 = nothing,
    y0 = nothing,
    nonvolatiles = nothing,
    f_limit = 0.0,
    atol = 1e-8,
    rtol = 1e-12,
    max_iters = 10^4)

    if T0 == y0 == vol0 == nothing
        return ChemPotBubbleTemperature{Nothing}(vol0,T0,y0,nonvolatiles,f_limit,atol,rtol,max_iters)
    elseif (T0 == y0 == nothing) && !isnothing(vol0)
        vl,vv = promote(vol0[1],vol0[2])
        return ChemPotBubbleTemperature{typeof(vl)}(vol0,T0,y0,nonvolatiles,f_limit,atol,rtol,max_iters)
    elseif (vol0 == y0 == nothing) && !isnothing(T0)
        T0 = float(T0)
        return ChemPotBubbleTemperature{typeof(T0)}(vol0,T0,y0,nonvolatiles,f_limit,atol,rtol,max_iters)
    elseif (T0 == vol0 == nothing) && !isnothing(y0)
        T = eltype(y0)
        return ChemPotBubbleTemperature{T}(vol0,T0,y0,nonvolatiles,f_limit,atol,rtol,max_iters)
    elseif !isnothing(vol0) && !isnothing(T0) && !isnothing(y0)
        vl,vv,T0,_ = promote(vol0[1],vol0[2],T0,first(y0))
        T = eltype(vl)
        y0 = convert(Vector{T},y0)
        return ChemPotBubbleTemperature{T}(vol0,T0,y0,nonvolatiles,f_limit,atol,rtol,max_iters)
    elseif !isnothing(vol0) && !isnothing(y0)
        vl,vv,_ = promote(vol0[1],vol0[2],first(y0))
        T = eltype(vl)
        y0 = convert(Vector{T},y0)
        return ChemPotBubbleTemperature{T}(vol0,T0,y0,nonvolatiles,f_limit,atol,rtol,max_iters)
    elseif  !isnothing(T0) && !isnothing(y0)
        T0,_ = promote(T0,first(y0))
        T = eltype(T0)
        y0 = convert(Vector{T},y0)
        return ChemPotBubbleTemperature{T}(vol0,T0,y0,nonvolatiles,f_limit,atol,rtol,max_iters)
    else
        throw(error("invalid specification for bubble temperature"))
    end
end

function bubble_temperature_impl(model::EoSModel,p,x,method::ChemPotBubbleTemperature)
    T0,vl,vv,y0 = bubble_temperature_init(model,p,x,method.vol0,method.T0,method.y0)
    if !isnothing(method.nonvolatiles)
        volatiles = [!in(x,method.nonvolatiles) for x in model.components]
        model_y,volatiles = index_reduction(model,volatiles)
        y0 = y0[volatiles]
        ts = T_scales(model_y,y0)
    else
        volatiles = fill(true,length(model))
        model_y = nothing
        ts = T_scales(model)
    end
    v0 = vcat(T0,log10(vl),log10(vv),y0[1:end-1])    
    pmix = p_scale(model,x)
    f!(F,z) = Obj_bubble_temperature(model,model_y, F, p, z[1], exp10(z[2]), exp10(z[3]), x, z[4:end],ts,pmix,volatiles)
    r  =Solvers.nlsolve(f!,v0,LineSearch(Newton()),NLSolvers.NEqOptions(method))
    sol = Solvers.x_sol(r)
    T   = sol[1]
    v_l = exp10(sol[2])
    v_v = exp10(sol[3])
    y_r = FractionVector(sol[4:end])
    y = index_expansion(y_r,volatiles)
    return T, v_l, v_v, y
end

function Obj_bubble_temperature(model::EoSModel,model_y, F, p, T, v_l, v_v, x, y,ts,ps,_view)
    F = μp_equality(model::EoSModel, model_y, F, T, v_l, v_v, x, FractionVector(y),ts,ps,_view)
    F[end] = (pressure(model,v_l,T,x) - p)/ps
    return F
end

#used by LLE_temperature
function Obj_bubble_temperature(model::EoSModel, F, p, T, v_l, v_v, x, y,ts,ps)
    return Obj_bubble_temperature(model,nothing, F, p, T, v_l, v_v, x, y,ts,ps,nothing)
end

export ChemPotBubblePressure, ChemPotBubbleTemperature