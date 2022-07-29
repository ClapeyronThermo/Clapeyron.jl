## Bubble pressure solver
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
        vl,vv,p0,_ = promote(vol0[1],vol0[2],p0,first(y))
        T = eltype(vl)
        y0 = convert(Vector{T},y0)
        return ChemPotBubblePressure{T}(vol0,p0,y0,nonvolatiles,f_limit,atol,rtol,max_iters)
    elseif !isnothing(vol0) && !isnothing(y0)
        vl,vv,_ = promote(vol0[1],vol0[2],first(y))
        T = eltype(vl)
        y0 = convert(Vector{T},y0)
        return ChemPotBubblePressure{T}(vol0,p0,y0,nonvolatiles,f_limit,atol,rtol,max_iters)
    elseif  !isnothing(p0) && !isnothing(y0)
        p0,_ = promote(p0,first(y))
        T = eltype(p0)
        y0 = convert(Vector{T},y0)
        return ChemPotBubblePressure{T}(vol0,p0,y0,nonvolatiles,f_limit,atol,rtol,max_iters)
    else
        throw(error("invalid specification for bubble pressure"))
    end
end

function bubble_pressure_impl(model::EoSModel, T, x,method::ChemPotBubblePressure)
    p0,vl,vv,y0 = bubble_pressure_init(model,T,x,method.vol0,method.p0,method.y0)
    v0 = vcat(log10(vl),log10(vv),y0[1:end-1])
    #xcache = zeros(eltype(x0),len)
    len = length(v0)
    Fcache = zeros(eltype(v0),len)
    model_r,x_r = model,x
    ts = T_scales(model_r,x_r)
    pmix = p_scale(model_r,x_r)
    f!(F,z) = Obj_bubble_pressure(model, F, T, exp10(z[1]),exp10(z[2]),x,z[3:end],ts,pmix)
    r  =Solvers.nlsolve(f!,v0,LineSearch(Newton()))
    sol = Solvers.x_sol(r)
    v_l = exp10(sol[1])
    v_v = exp10(sol[2])
    y = FractionVector(sol[3:end])
    P_sat = pressure(model,v_l,T,x_r)
    return (P_sat, v_l, v_v, y)
end


function Obj_bubble_pressure(model::EoSModel, F, T, v_l, v_v, x, y,ts,ps)
    return μp_equality(model::EoSModel, F, T, v_l, v_v, x, FractionVector(y),ts,ps)
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
        vl,vv,T0,_ = promote(vol0[1],vol0[2],T0,first(y))
        T = eltype(vl)
        y0 = convert(Vector{T},y0)
        return ChemPotBubbleTemperature{T}(vol0,T0,y0,nonvolatiles,f_limit,atol,rtol,max_iters)
    elseif !isnothing(vol0) && !isnothing(y0)
        vl,vv,_ = promote(vol0[1],vol0[2],first(y))
        T = eltype(vl)
        y0 = convert(Vector{T},y0)
        return ChemPotBubbleTemperature{T}(vol0,T0,y0,nonvolatiles,f_limit,atol,rtol,max_iters)
    elseif  !isnothing(T0) && !isnothing(y0)
        T0,_ = promote(T0,first(y))
        T = eltype(T0)
        y0 = convert(Vector{T},y0)
        return ChemPotBubbleTemperature{T}(vol0,T0,y0,nonvolatiles,f_limit,atol,rtol,max_iters)
    else
        throw(error("invalid specification for bubble temperature"))
    end
end

function bubble_temperature_impl(model::EoSModel,p,x,method::ChemPotBubbleTemperature)
    model_r,x_r = model,x
    #x_r = x[idx_r]
    ts = T_scales(model_r)
    pmix = p_scale(model_r,x_r)
    T0,vl,vv,y0 = bubble_temperature_init(model,p,x,method.vol0,method.T0,method.y0)
    v0 = vcat(T0,log10(vl),log10(vv),y0[1:end-1])
    len = length(v0)
    f!(F,z) = Obj_bubble_temperature(model_r, F, p, z[1], exp10(z[2]), exp10(z[3]), x_r, z[4:end],ts,pmix)
    r  =Solvers.nlsolve(f!,v0,LineSearch(Newton()))
    sol = Solvers.x_sol(r)
    T   = sol[1]
    v_l = exp10(sol[2])
    v_v = exp10(sol[3])
    y_r = FractionVector(sol[4:end])
    y = zeros(length(model))
    y = y_r
    #y[idx_r] = y_r
    return T, v_l, v_v, y
end

function Obj_bubble_temperature(model::EoSModel, F, p, T, v_l, v_v, x, y,ts,ps)
    F = μp_equality(model::EoSModel, F, T, v_l, v_v, x, FractionVector(y),ts,ps)
    F[end] = (pressure(model,v_l,T,x) - p)/ps
    return F
end