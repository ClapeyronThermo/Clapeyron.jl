struct ActivityQT{T} <: FlashMethod
    data::FugEnum.BubbleDew
    vol0::Union{Nothing,Tuple{T,T}}
    p0::Union{Nothing,T}
    w0::Union{Nothing,Vector{T}}
    non_in_w::Union{Nothing,Vector{String}}
    itmax_ss::Int64
    rtol_ss::Float64
    verbose::Bool
end

function Solvers.primalval(method::ActivityQT{T}) where T
    if T == Nothing
        return Solvers.primalval_struct(method,T)
    else
        return Solvers.primalval_struct(method,Solvers.primal_eltype(T))
    end
end

function index_reduction(method::ActivityQT,idx_r)
    if !isnothing(method.w0)
        method_r = deepcopy(method)
        w0_new = method.w0[idx_r]
        resize!(method_r.w0,length(w0_new))
        method_r.w0 .= w0_new
        return method_r
    end
    return method
end

function bdt_flash_impl(wrapper::PTFlashWrapper,T,z,method::ActivityQT)
    model = wrapper.model
    TT = Base.promote_eltype(wrapper,T,x)
    in_equilibria = comps_in_equilibria(component_list(wrapper),method.non_in_w)

    if is_lle(method.equilibrium)
        phasex,phasey = :liquid,:liquid
        p,volx,voly,y = LLE_pressure_init(wrapper,T,z,method.vol0,method.p0,method.w0,in_equilibria)
        x = similar(y)
        x .= z
    elseif bubble
        phasex,phasey = :liquid,:vapour
        p,volx,voly,y = bubble_pressure_init(wrapper,T,z,method.vol0,method.p0,method.w0,in_equilibria)
        x = similar(y)
        x .= z
    else
        phasex,phasey = :liquid,:vapour
        p,volx,voly,x = dew_pressure_init(wrapper,T,z,method.vol0,method.p0,method.w0,in_equilibria)
        y = similar(x)
        y .= z
    end

    cache = ∂lnϕ_cache(model, p, T, x, Val{false}())
    lnK = similar(x,TT)
    pold = -p
    lnK = similar(x,TT)
    K = similar(x,TT)
    for k in 1:method.itmax_ss
        lnϕx, volx = modified_lnϕ(wrapper, p, T, x, cache; phase = phasex, vol0 = volx)
        lnK .= lnϕx
        lnϕy, voly = modified_lnϕ(wrapper, p, T, y, cache; phase = phasey, vol0 = voly)
        lnK .-= lnϕy
        lnK .+= log(p)
        K .= exp.(lnK)
        if bubble
            for i in 1:length(K)
                in_equilibria[i] || (K[i] = 0)
            end
            y .= K .* x
            pold = p
            p = sum(y)
            y ./= p
        else
            for i in 1:length(K)
                in_equilibria[i] || (K[i] = Inf)
            end
            x .= y ./ K
            pold = p
            p = 1/sum(x)
            x .*= p
        end
        err = abs(pold-p)/p
        !isfinite(err) && break
        err < method.rtol_ss && break
    end
    if iszero(volx)
        volx = volume(wrapper,p,T,x,phase = phasex)
    end

    if iszero(voly)
        voly = volume(wrapper,p,T,y,phase = phasey)
    end

    return (p,volx,voly,y)
end

function ActivityQT(data::FugEnum.BubbleDew;vol0 = nothing,
                                p0 = nothing,
                                w0 = nothing,
                                non_in_w = nothing,
                                itmax_ss = 40,
                                rtol_ss = 1e-8,
                                verbose = false)

    if p0 == w0 == vol0 == nothing
        return ActivityQT{Nothing}(data,vol0,p0,w0,non_in_w,itmax_ss,rtol_ss,verbose)
    elseif (p0 == w0 == nothing) && !isnothing(vol0)
        vl,vv = promote(vol0[1],vol0[2])
        return ActivityQT{typeof(vl)}(data,vol0,p0,w0,non_in_w,itmax_ss,rtol_ss,verbose)
    elseif (vol0 == w0 == nothing) && !isnothing(p0)
        p0 = float(p0)
        return ActivityQT{typeof(p0)}(data,vol0,p0,w0,non_in_w,itmax_ss,rtol_ss,verbose)
    elseif (p0 == vol0 == nothing) && !isnothing(w0)
        T = eltype(w0)
        return ActivityQT{T}(data,vol0,p0,w0,non_in_w,itmax_ss,rtol_ss,verbose)
    elseif !isnothing(vol0) && !isnothing(p0) && !isnothing(w0)
        vl,vv,p0,_ = promote(vol0[1],vol0[2],p0,first(w0))
        T = eltype(vl)
        w0 = convert(Vector{T},w0)
        return ActivityQT{T}(data,vol0,p0,w0,non_in_w,itmax_ss,rtol_ss,verbose)
    elseif !isnothing(vol0) && !isnothing(w0)
        vl,vv,_ = promote(vol0[1],vol0[2],first(w0))
        T = eltype(vl)
        w0 = convert(Vector{T},w0)
        return ActivityQT{T}(data,vol0,p0,w0,non_in_w,itmax_ss,rtol_ss,verbose)
    elseif  !isnothing(p0) && !isnothing(w0)
        p0,_ = promote(p0,first(w0))
        T = eltype(p0)
        w0 = convert(Vector{T},w0)
        return ActivityQT{T}(data,vol0,p0,w0,non_in_w,itmax_ss,rtol_ss,verbose)
    else
        invalid_bd_error(data)
    end
end


## Bubble pressure solver
"""
    ActivityBubblePressure(kwargs...)

Function to compute [`bubble_pressure`](@ref) using Activity Coefficients.
On activity coefficient models it solves the problem via succesive substitucion.
On Helmholtz-based models, it approximates the activity coefficient using the saturated pure state as reference.

Inputs:
- `y0 = nothing`: optional, initial guess for the vapor phase composition
- `p0 = nothing`: optional, initial guess for the bubble pressure `[Pa]`
- `vol0 = nothing`: optional, initial guesses for the liquid and vapor phase volumes `[m³]`
- `atol = 1e-8`: optional, absolute tolerance of the non linear system of equations
- `rtol = 1e-12`: optional, relative tolerance of the non linear system of equations
- `itmax_ss = 40`: optional, maximum number of sucesive substitution iterations
- `nonvolatiles = nothing`: optional, Vector of strings containing non volatile compounds. Those will be set to zero on the vapour phase.
"""
function ActivityBubblePressure(;vol0 = nothing,
                                p0 = nothing,
                                y0 = nothing,
                                nonvolatiles = nothing,
                                itmax_ss = 40,
                                rtol_ss = 1e-8,
                                verbose = false)

    w0 = y0
    non_in_w = nonvolatiles
    return ActivityQT(FugEnum.BUBBLE_PRESSURE,vol0,p0,w0,non_in_w,itmax_ss,rtol_ss,verbose)
end

"""
    ActivityDewPressure(kwargs...)

Function to compute [`dew_pressure`](@ref) using Activity Coefficients.
On activity coefficient models it solves the problem via succesive substitucion.
On Helmholtz-based models, it approximates the activity coefficient using the saturated pure state as reference.

Inputs:
- `x0 = nothing`: optional, initial guess for the liquid phase composition
- `p0 = nothing`: optional, initial guess for the dew pressure `[Pa]`
- `vol0 = nothing`: optional, initial guesses for the liquid and vapor phase volumes `[m³]`
- `atol = 1e-8`: optional, absolute tolerance of the non linear system of equations
- `rtol = 1e-12`: optional, relative tolerance of the non linear system of equations
- `itmax_ss = 40`: optional, maximum number of sucesive substitution iterations
- `noncondensables`: optional, Vector of strings containing non condensable compounds. Those will be set to zero on the liquid phase.

"""
function ActivityDewPressure(;vol0 = nothing,
                                p0 = nothing,
                                x0 = nothing,
                                noncondensables = nothing,
                                itmax_ss = 40,
                                rtol_ss = 1e-8,
                                verbose = false)

    w0 = x0
    non_in_w = noncondensables
    return ActivityQT(FugEnum.DEW_PRESSURE,vol0,p0,w0,non_in_w,itmax_ss,rtol_ss,verbose)
end
