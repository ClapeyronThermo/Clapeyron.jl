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

"""
    bubble_temperature(model::EoSModel, p, x; T0 = x0_bubble_pressure(model,p,x))

calculates the bubble temperature and properties at a given pressure.
Returns a tuple, containing:
- Bubble Temperature `[K]`
- liquid volume at Bubble Point [`m³`]
- vapour volume at Bubble Point [`m³`]
- Gas composition at Bubble Point
"""
function bubble_temperature(model::EoSModel,p,x;v0=nothing)
    model_r,idx_r = index_reduction(model,x)
    if length(model_r)==1
        (T_sat,v_l,v_v) = saturation_temperature(model_r,p,v0)
        return (T_sat,v_l,v_v,x)
    end
    x_r = x[idx_r]
    ts = T_scales(model_r)
    pmix = p_scale(model_r,x_r)
    if v0 === nothing
        v0 = x0_bubble_temperature(model_r,p,x_r)
    end
    
    len = length(v0[1:end-1])
    Fcache = zeros(eltype(v0[1:end-1]),len)
    f!(F,z) = Obj_bubble_temperature(model_r, F, p, z[1], exp10(z[2]), exp10(z[3]), x_r, z[4:end],ts,pmix)
    r  =Solvers.nlsolve(f!,v0[1:end-1],LineSearch(Newton()))
    sol = Solvers.x_sol(r)
    T   = sol[1]
    v_l = exp10(sol[2])
    v_v = exp10(sol[3])
    y_r = FractionVector(sol[4:end])
    y = zeros(length(model))
    y[idx_r] = y_r
    return T, v_l, v_v, y
end

function Obj_bubble_temperature(model::EoSModel, F, p, T, v_l, v_v, x, y,ts,ps)
    F = μp_equality(model::EoSModel, F, T, v_l, v_v, x, FractionVector(y),ts,ps)
    F[end] = (pressure(model,v_l,T,x) - p)/ps
    return F
end

function x0_bubble_temperature(model::EoSModel,p,x)
    comps = length(model)   
    pure = split_model(model)
    crit = crit_pure.(pure)
    
    p_c = [tup[2] for tup in crit]
    V_c = [tup[3] for tup in crit]
    _0 = zero(p+first(x))
    replaceP = p_c .< p
    T_sat = fill(_0,comps)
    V_l_sat = fill(_0,comps)
    V_v_sat = fill(_0,comps)
    for i in 1:comps
        crit_i = crit[i]
        Tci,Pci,Vci = crit_i
        if !replaceP[i]
            Ti,Vli,Vvi = saturation_temperature(pure[i],p,AntoineSaturation(crit = crit_i))
        else
            
            Ti,Vli,Vvi = Tci,Vci,1.2*Vci  
        end
        T_sat[i] = Ti
        V_l_sat[i] = Vli
        V_v_sat[i] = Vvi
    end    
    Tb = extrema(T_sat).*(0.9,1.1)

    V0_l = zero(p)
    V0_v = zero(p)
    f(T) = antoine_bubble(pure,T,x,crit)[1]-p
    fT = Roots.ZeroProblem(f,Tb)

    T0 = Roots.solve(fT,Roots.Order0())
    p,y = antoine_bubble(pure,T0,x,crit)
    for i in 1:length(x)
        if !replaceP[i]
            V0_v += y[i]*V_v_sat[i]
            V0_l += x[i]*V_l_sat[i]
        else 
            V0_v += y[i]*V_c[i]*1.2
            V0_l += x[i]*V_c[i]
        end
    end
    prepend!(y,log10.([V0_l,V0_v]))
    prepend!(y,T0)

    return y
end

function antoine_bubble(pure,T,x,crit)
    pᵢ = aprox_psat.(pure,T,crit)
    p = sum(x.*pᵢ)
    y = x.*pᵢ./p
    ysum = 1/∑(y)
    y    = y.*ysum
    return p,y
end

