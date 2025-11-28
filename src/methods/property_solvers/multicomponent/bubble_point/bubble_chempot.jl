## Bubble pressure solver
"""
    ChemPotBubblePressure(kwargs...)

Function to compute [`bubble_pressure`](@ref) via chemical potentials.
It directly solves the equality of chemical potentials system of equations.

Inputs:
- `y0 = nothing`: optional, initial guess for the vapor phase composition
- `p0 = nothing`: optional, initial guess for the bubble pressure `[Pa]`
- `vol0 = nothing`: optional, initial guesses for the liquid and vapor phase volumes `[m³]`
- `atol = 1e-8`: optional, absolute tolerance of the non linear system of equations
- `rtol = 1e-12`: optional, relative tolerance of the non linear system of equations
- `max_iters = 1000`: optional, maximum number of iterations
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
    ss::Bool
end

function ChemPotBubblePressure(;vol0 = nothing,
                                p0 = nothing,
                                y0 = nothing,
                                nonvolatiles = nothing,
                                f_limit = 0.0,
                                atol = 1e-8,
                                rtol = 1e-12,
                                max_iters = 1000,
                                ss = false)

    if p0 == y0 == vol0 == nothing
        return ChemPotBubblePressure{Nothing}(vol0,p0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,ss)
    elseif (p0 == y0 == nothing) && !isnothing(vol0)
        vl,vv = promote(vol0[1],vol0[2])
        return ChemPotBubblePressure{typeof(vl)}(vol0,p0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,ss)
    elseif (vol0 == y0 == nothing) && !isnothing(p0)
        p0 = float(p0)
        return ChemPotBubblePressure{typeof(p0)}(vol0,p0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,ss)
    elseif (p0 == vol0 == nothing) && !isnothing(y0)
        T = eltype(y0)
        return ChemPotBubblePressure{T}(vol0,p0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,ss)
    elseif !isnothing(vol0) && !isnothing(p0) && !isnothing(y0)
        vl,vv,p0,_ = promote(vol0[1],vol0[2],p0,first(y0))
        T = eltype(vl)
        y0 = convert(Vector{T},y0)
        return ChemPotBubblePressure{T}(vol0,p0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,ss)
    elseif !isnothing(vol0) && !isnothing(y0)
        vl,vv,_ = promote(vol0[1],vol0[2],first(y0))
        T = eltype(vl)
        y0 = convert(Vector{T},y0)
        return ChemPotBubblePressure{T}(vol0,p0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,ss)
    elseif  !isnothing(p0) && !isnothing(y0)
        p0,_ = promote(p0,first(y0))
        T = eltype(p0)
        y0 = convert(Vector{T},y0)
        return ChemPotBubblePressure{T}(vol0,p0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,ss)
    else
        throw(error("invalid specification for bubble pressure"))
    end
end

function bubble_pressure_impl(model::EoSModel, T, x,method::ChemPotBubblePressure)
    volatiles = comps_in_equilibria(component_list(model),method.nonvolatiles)
    p0,vl0,vv0,y00 = bubble_pressure_init(model,T,x,method.vol0,method.p0,method.y0,volatiles)
    is_non_volatile = !isnothing(method.nonvolatiles)
    model_y,_ = index_reduction(model,volatiles)
    y0 = y00[volatiles]
    data = FugEnum.BUBBLE_PRESSURE
    
    neq = count(volatiles)
    cache = Clapeyron.fug_bubbledew_cache(model,model_y,T,T,x,x,Val{false}())
    y_r,_,_,_,_,_,p_cache,_,_ = cache
    inc0 = similar(y_r,length(y_r)+2)
    inc0[1:neq] .= log.(y0 ./ @view(x[volatiles]))
    inc0[neq+1] = log(vl0)
    inc0[neq+2] = log(vv0)
    opts = NLSolvers.NEqOptions(method)
    if true #is_non_volatile
        problem = _mu_OF_neq(model,T,x,data,cache)
        sol = Solvers.nlsolve(problem, inc0, Solvers.LineSearch(Solvers.Newton2(inc0),Static(1.0)),opts)
        inc = Solvers.x_sol(sol)
        !all(<(sol.options.f_abstol),sol.info.best_residual) && (inc .= NaN)
    else
        problem = _mu_OF_neq(model,model_y,T,x,volatiles,data,cache)
        sol_vol = Solvers.nlsolve(problem, inc0, Solvers.LineSearch(Solvers.Newton2(inc0),Static(1.0)),opts)
        inc = Solvers.x_sol(sol_vol)
        !all(<(sol_vol.options.f_abstol),sol_vol.info.best_residual) && (inc .= NaN)
    end
    lnK = @view inc[1:end-2]
    y_r .= exp.(lnK) .* @view(x[volatiles])
    y = index_expansion(y_r,volatiles)
    vl = exp(inc[end-1])
    vv = exp(inc[end])
    px,py = p_cache[]
    p = 0.5*(px + py)
    return (p, vl, vv, y)
end


function bubble_pressure2(model::EoSModel, T, x,method::ChemPotBubblePressure)
    volatiles = comps_in_equilibria(component_list(model),method.nonvolatiles)
    p0,vl,vv,y0 = bubble_pressure_init(model,T,x,method.vol0,method.p0,method.y0,volatiles)
    is_non_volatile = !isnothing(method.nonvolatiles)
    model_y,_ = index_reduction(model,volatiles)
    y0 = y0[volatiles]

    ηl = η_from_v(model, vl, T, x)
    if is_non_volatile
        ηv = η_from_v(model_y, vv, T, y0)
    else
        ηv = η_from_v(model, vv, T, y0)
    end

    # select component with highest fraction as pivot
    idx_max = argmax(y0)
    v0 = Vector{eltype(y0)}(undef, 2 + length(y0) - 1)
    v0[1],v0[2] = ηl,ηv
    copy_without_pivot!(view(v0, 3:lastindex(v0)), y0, idx_max)
    f! = let model = model, model_y = model_y, T=T, x=x, volatiles=volatiles, idx_max=idx_max
        (F,z) -> Obj_bubble_pressure(model, model_y, F, T, z[1], z[2], x, @view(z[3:end]), volatiles, idx_max)
    end
    r = Solvers.nlsolve(f!, v0,
        LineSearch(Newton2(v0)),
        NLSolvers.NEqOptions(method),
        ForwardDiff.Chunk{min(length(v0), 8)}()
    )
    sol = Solvers.x_sol(r)
    !all(<(r.options.f_abstol),r.info.best_residual) && (sol .= NaN)
    v_l = v_from_η(model,sol[1],T,x)
    y_r = FractionVector(@view(sol[3:end]),idx_max)
    v_v = v_from_η(model,model_y,sol[2],T,y_r)
    y_sol = index_expansion(y_r,volatiles)
    P_sat = pressure(model,v_l,T,x)
    return (P_sat, v_l, v_v, y_sol)
end


function Obj_bubble_pressure(model::EoSModel, model_y, F, T, ηl, ηv, x, y, _view,yy_i)
    v_l = v_from_η(model,ηl,T,x)
    yy = FractionVector(y,yy_i)
    v_v = v_from_η(model_y,ηv,T,yy)
    v = (v_l,v_v)
    w = (x,yy)
    if isnothing(_view) || all(_view)
        return μp_equality2(model, nothing, F, Tspec(T), v, w, _view)
    else
        return μp_equality2(model, model_y, F, Tspec(T), v, w, _view)
    end
end

#used by LLE_pressure
function Obj_bubble_pressure(model::EoSModel, F, T, ηl, ηv, x, y)
    return Obj_bubble_pressure(model, model, F, T, ηl, ηv, x, y, nothing, length(model))
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
    ss::Bool
end

"""
    ChemPotBubbleTemperature(kwargs...)

Function to compute [`bubble_temperature`](@ref) via chemical potentials.
It directly solves the equality of chemical potentials system of equations.

Inputs:
- `y0 = nothing`: optional, initial guess for the vapor phase composition.
- `T0 = nothing`: optional, initial guess for the bubble temperature `[K]`.
- `vol0 = nothing`: optional, initial guesses for the liquid and vapor phase volumes `[m³]`
- `atol = 1e-8`: optional, absolute tolerance of the non linear system of equations
- `rtol = 1e-12`: optional, relative tolerance of the non linear system of equations
- `max_iters = 1000`: optional, maximum number of iterations
- `nonvolatiles = nothing`: optional, Vector of strings containing non volatile compounds. those will be set to zero on the vapour phase.
"""
function ChemPotBubbleTemperature(;vol0 = nothing,
    T0 = nothing,
    y0 = nothing,
    nonvolatiles = nothing,
    f_limit = 0.0,
    atol = 1e-8,
    rtol = 1e-12,
    max_iters = 10^3,
    ss = false)

    if T0 == y0 == vol0 == nothing
        return ChemPotBubbleTemperature{Nothing}(vol0,T0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,ss)
    elseif (T0 == y0 == nothing) && !isnothing(vol0)
        vl,vv = promote(vol0[1],vol0[2])
        return ChemPotBubbleTemperature{typeof(vl)}(vol0,T0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,ss)
    elseif (vol0 == y0 == nothing) && !isnothing(T0)
        T0 = float(T0)
        return ChemPotBubbleTemperature{typeof(T0)}(vol0,T0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,ss)
    elseif (T0 == vol0 == nothing) && !isnothing(y0)
        T = eltype(y0)
        return ChemPotBubbleTemperature{T}(vol0,T0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,ss)
    elseif !isnothing(vol0) && !isnothing(T0) && !isnothing(y0)
        vl,vv,T0,_ = promote(vol0[1],vol0[2],T0,first(y0))
        T = eltype(vl)
        y0 = convert(Vector{T},y0)
        return ChemPotBubbleTemperature{T}(vol0,T0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,ss)
    elseif !isnothing(vol0) && !isnothing(y0)
        vl,vv,_ = promote(vol0[1],vol0[2],first(y0))
        T = eltype(vl)
        y0 = convert(Vector{T},y0)
        return ChemPotBubbleTemperature{T}(vol0,T0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,ss)
    elseif  !isnothing(T0) && !isnothing(y0)
        T0,_ = promote(T0,first(y0))
        T = eltype(T0)
        y0 = convert(Vector{T},y0)
        return ChemPotBubbleTemperature{T}(vol0,T0,y0,nonvolatiles,f_limit,atol,rtol,max_iters,ss)
    else
        throw(error("invalid specification for bubble temperature"))
    end
end

function bubble_temperature_impl(model::EoSModel,p,x,method::ChemPotBubbleTemperature)
    is_non_volatile = !isnothing(method.nonvolatiles)
    volatiles = comps_in_equilibria(component_list(model),method.nonvolatiles)
    model_y,_ = index_reduction(model,volatiles)
    T0,vl,vv,y0 = bubble_temperature_init(model,p,x,method.vol0,method.T0,method.y0,volatiles)
    y0 = y0[volatiles]

    ηl = η_from_v(model, vl, T0, x)
    if is_non_volatile
        ηv = η_from_v(model_y, vv, T0, y0)
    else
        ηv = η_from_v(model, vv, T0, y0)
    end

    # select component with highest fraction as pivot
    idx_max = argmax(y0)
    v0 = Vector{eltype(y0)}(undef, 3 + length(y0) - 1)
    v0[1],v0[2],v0[3] = T0,ηl,ηv
    copy_without_pivot!(view(v0, 4:lastindex(v0)), y0, idx_max)
    f! = let model = model, model_y = model_y, p=p, x=x, volatiles=volatiles, idx_max=idx_max
        (F,z) -> Obj_bubble_temperature(model, model_y, F, p, z[1], z[2], z[3], x, @view(z[4:end]), volatiles, idx_max)
    end
    r = Solvers.nlsolve(f!, v0,
        LineSearch(Newton2(v0)), 
        NLSolvers.NEqOptions(method),
        ForwardDiff.Chunk{min(length(v0), 8)}()
    )
    sol = Solvers.x_sol(r)
    !all(<(r.options.f_abstol),r.info.best_residual) && (sol .= NaN)
    T = sol[1]
    y_r = FractionVector(@view(sol[4:end]), idx_max)
    v_l = v_from_η(model, sol[2], T, x)
    v_v = v_from_η(model, model_y, sol[3], T, y_r)
    y = index_expansion(y_r,volatiles)
    return T, v_l, v_v, y
end

function Obj_bubble_temperature(model::EoSModel, model_y, F, p, T, ηl, ηv, x, y, _view, yy_i)
    yy = FractionVector(y,yy_i)
    vl = v_from_η(model, ηl, T, x)
    vv = v_from_η(model_y, ηv, T, yy)
    v = (vl,vv)
    w = (x,yy)
    if isnothing(_view) || all(_view)
        return μp_equality2(model, nothing, F, Pspec(p,T), v, w, _view)
    else
        return μp_equality2(model, model_y, F, Pspec(p,T), v, w, _view)
    end
end

#used by LLE_temperature
function Obj_bubble_temperature(model::EoSModel, F, p, T, ηl, ηv, x, y)
    return Obj_bubble_temperature(model,model, F, p, T, ηl, ηv, x, y,nothing,length(model))
end

export ChemPotBubblePressure, ChemPotBubbleTemperature