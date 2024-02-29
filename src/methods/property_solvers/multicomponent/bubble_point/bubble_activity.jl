## Bubble pressure solver
"""
    ActivityBubblePressure(kwargs...)

Function to compute [`bubble_pressure`](@ref) using Activity Coefficients.
On activity coefficient models it solves the problem via succesive substitucion.
On helmholtz-based models, it uses the Chapman approximation for activity coefficients.

Inputs:
- `gas_fug = true`: if the solver uses gas fugacity coefficients. on `ActivityModel` is set by default to `false`
- `poynting = true`: if the solver use the poynting correction on the liquid fugacity coefficients. on `ActivityModel` is set by default to `false`
- `y0 = nothing`: optional, initial guess for the vapor phase composition
- `p0 = nothing`: optional, initial guess for the bubble pressure [`Pa`]
- `vol0 = nothing`: optional, initial guesses for the liquid and vapor phase volumes
- `atol = 1e-8`: optional, absolute tolerance of the non linear system of equations
- `rtol = 1e-12`: optional, relative tolerance of the non linear system of equations
- `itmax_ss = 40`: optional, maximum number of sucesive substitution iterations
"""
struct ActivityBubblePressure{T} <: BubblePointMethod
    vol0::Union{Nothing,Tuple{T,T}}
    p0::Union{Nothing,T}
    y0::Union{Nothing,Vector{T}}
    nonvolatiles::Union{Nothing,Vector{String}}
    itmax_ss::Int64
    rtol_ss::Float64
    gas_fug::Bool
    poynting::Bool
end

function ActivityBubblePressure(;vol0 = nothing,
                                p0 = nothing,
                                y0 = nothing,
                                nonvolatiles = nothing,
                                itmax_ss = 40,
                                rtol_ss = 1e-8,
                                gas_fug = true,
                                poynting = true)

    if p0 == y0 == vol0 == nothing
        return ActivityBubblePressure{Nothing}(vol0,p0,y0,nonvolatiles,itmax_ss,rtol_ss,gas_fug,poynting)
    elseif (p0 == y0 == nothing) && !isnothing(vol0)
        vl,vv = promote(vol0[1],vol0[2])
        return ActivityBubblePressure{typeof(vl)}(vol0,p0,y0,nonvolatiles,itmax_ss,rtol_ss,gas_fug,poynting)
    elseif (vol0 == y0 == nothing) && !isnothing(p0)
        p0 = float(p0)
        return ActivityBubblePressure{typeof(p0)}(vol0,p0,y0,nonvolatiles,itmax_ss,rtol_ss,gas_fug,poynting)
    elseif (p0 == vol0 == nothing) && !isnothing(y0)
        T = eltype(y0)
        return ActivityBubblePressure{T}(vol0,p0,y0,nonvolatiles,itmax_ss,rtol_ss,gas_fug,poynting)
    elseif !isnothing(vol0) && !isnothing(p0) && !isnothing(y0)
        vl,vv,p0,_ = promote(vol0[1],vol0[2],p0,first(y0))
        T = eltype(vl)
        y0 = convert(Vector{T},y0)
        return ActivityBubblePressure{T}(vol0,p0,y0,nonvolatiles,itmax_ss,rtol_ss,gas_fug,poynting)
    elseif !isnothing(vol0) && !isnothing(y0)
        vl,vv,_ = promote(vol0[1],vol0[2],first(y0))
        T = eltype(vl)
        y0 = convert(Vector{T},y0)
        return ActivityBubblePressure{T}(vol0,p0,y0,nonvolatiles,itmax_ss,rtol_ss,gas_fug,poynting)
    elseif  !isnothing(p0) && !isnothing(y0)
        p0,_ = promote(p0,first(y0))
        T = eltype(p0)
        y0 = convert(Vector{T},y0)
        return ActivityBubblePressure{T}(vol0,p0,y0,nonvolatiles,itmax_ss,rtol_ss,gas_fug,poynting)
    else
        throw(error("invalid specification for bubble pressure"))
    end
end


function bubble_pressure_impl(model,T,x,method::ActivityBubblePressure)
    R̄ = Rgas(model)
    pure = split_model(model)
    sat = saturation_pressure.(pure,T)
    p_pure = first.(sat)
    vl_pure = getindex.(sat,2)
    vv_pure = last.(sat)

    _0  =zero(eltype(p_pure))
    nan = _0/_0

    if isnothing(method.vol0)
        vl = dot(vl_pure,x)
        vv = zero(vl)
    else
        vl,vv = method.vol0
    end
   
    if isnan(vl)
        return _0,_0,_0,x
    end

    if any(isnan,p_pure)
        _0  =zero(eltype(p_pure))
        nan = 
        return _0,_0,_0,x
    end

    if isnothing(method.p0)
        pmix = pressure(model,vl,T,x)
    else
        pmix = method.p0
    end

    μmix = VT_chemical_potential_res(model,vl,T,x)
    ϕ = copy(μmix)
    y = copy(μmix)
    ϕpure = copy(μmix)
    ϕ .= 1
    ϕpure .= 1
    #y = x .* p_pure #raoult initialization
    #y ./= sum(y)
    #@show y
    RT = (R̄*T)
    if !isnothing(method.y0)
        y .= method.y0
        vv = dot(y,vv_pure)
        #if method.gas_fug
        #    μv = VT_chemical_potential_res!(ϕ,model,vv,T,y)
        #    ϕ .= exp.(μv ./ RT .- log.(pmix .* vv ./ RT))
        #end
    end


    pold = zero(pmix)
    γ = zeros(typeof(pmix),length(pure))
    #pure part
    μpure = only.(VT_chemical_potential_res.(pure,vl_pure,T))
    if method.gas_fug
        ϕpure .= exp.(μpure ./ RT .- log.(p_pure .* vl_pure ./ RT))
    end
    if method.poynting
        κ = VT_isothermal_compressibility.(pure,vl_pure,T)
    else
        κ = copy(ϕ)
        κ .= 0.0
    end
    for k in 1:method.itmax_ss
        for i in eachindex(γ)
            pᵢ = p_pure[i]
            vpureᵢ = vl_pure[i]
            μᵢ = μpure[i]
            ϕ̂ᵢ = ϕpure[i]
            γ[i] = exp(log(vpureᵢ/vl) + (μmix[i] - μᵢ)/RT -  vpureᵢ*(pmix -pᵢ)/RT)
            if method.poynting
                ln𝒫 = vpureᵢ*expm1(κ[i]*(pmix-pᵢ))/(κ[i]*RT) #see end of file
                𝒫 = exp(ln𝒫)
            else
                𝒫 = one(pᵢ)
            end
            #y[i]*ϕ[i]*P = x[i]*γ[i]*pᵢ*ϕ̂ᵢ*𝒫
            y[i] = x[i]*γ[i]*pᵢ*𝒫*ϕ̂ᵢ/ϕ[i] #really yᵢ*P, we normalize later
        end

        pold = pmix
        pmix = sum(y)
        y ./= pmix
        if iszero(vv)
            vv = dot(y,vv_pure)
        end
        vl = volume(model,pmix,T,x,vol0 = vl)
        if method.gas_fug
            logϕ, vv = lnϕ(model,pmix,T,y,phase = :vapor, vol0 = vv)
            ϕ .= exp.(logϕ)
        else
            vv = volume(model,pmix,T,y,phase =:vapor,vol0 = vv)
        end
        err = abs(pold-pmix)/pmix
        μmix = VT_chemical_potential_res!(μmix,model,vl,T,x)
        !isfinite(err) && break
        err < method.rtol_ss && break
    end
    return pmix,vl,vv,y
end

## Bubble pressure solver
"""
    ActivityBubbleTemperature(kwargs...)

Function to compute [`bubble_temperature`](@ref) using Activity Coefficients.
On activity coefficient models it solves the problem via succesive substitucion.
On helmholtz-based models, it uses the Chapman approximation for activity coefficients.

Inputs:
- `gas_fug = true`: if the solver uses gas fugacity coefficients. on `ActivityModel` is set by default to `false`
- `poynting = true`: if the solver use the poynting correction on the liquid fugacity coefficients. on `ActivityModel` is set by default to `false`
- `y0 = nothing`: optional, initial guess for the vapor phase composition
- `T0 = nothing`: optional, initial guess for the bubble temperature [`K`]
- `vol0 = nothing`: optional, initial guesses for the liquid and vapor phase volumes
- `atol = 1e-8`: optional, absolute tolerance of the non linear system of equations
- `rtol = 1e-12`: optional, relative tolerance of the non linear system of equations
- `itmax_ss = 40`: optional, maximum number of sucesive substitution iterations
"""
struct ActivityBubbleTemperature{T} <: BubblePointMethod
    vol0::Union{Nothing,Tuple{T,T}}
    T0::Union{Nothing,T}
    y0::Union{Nothing,Vector{T}}
    nonvolatiles::Union{Nothing,Vector{String}}
    itmax_ss::Int64
    rtol_ss::Float64
    gas_fug::Bool
    poynting::Bool
end

function ActivityBubbleTemperature(;vol0 = nothing,
                                T0 = nothing,
                                y0 = nothing,
                                nonvolatiles = nothing,
                                itmax_ss = 40,
                                rtol_ss = 1e-8,
                                gas_fug = true,
                                poynting = true)

    if T0 == y0 == vol0 == nothing
        return ActivityBubbleTemperature{Nothing}(vol0,T0,y0,nonvolatiles,itmax_ss,rtol_ss,gas_fug,poynting)
    elseif (T0 == y0 == nothing) && !isnothing(vol0)
        vl,vv = promote(vol0[1],vol0[2])
        return ActivityBubbleTemperature{typeof(vl)}(vol0,T0,y0,nonvolatiles,itmax_ss,rtol_ss,gas_fug,poynting)
    elseif (vol0 == y0 == nothing) && !isnothing(T0)
        T0 = float(T0)
        return ActivityBubbleTemperature{typeof(T0)}(vol0,T0,y0,nonvolatiles,itmax_ss,rtol_ss,gas_fug,poynting)
    elseif (T0 == vol0 == nothing) && !isnothing(y0)
        T = eltype(y0)
        return ActivityBubbleTemperature{T}(vol0,T0,y0,nonvolatiles,itmax_ss,rtol_ss,gas_fug,poynting)
    elseif !isnothing(vol0) && !isnothing(T0) && !isnothing(y0)
        vl,vv,T0,_ = promote(vol0[1],vol0[2],T0,first(y0))
        T = eltype(vl)
        y0 = convert(Vector{T},y0)
        return ActivityBubbleTemperature{T}(vol0,T0,y0,nonvolatiles,itmax_ss,rtol_ss,gas_fug,poynting)
    elseif !isnothing(vol0) && !isnothing(y0)
        vl,vv,_ = promote(vol0[1],vol0[2],first(y0))
        T = eltype(vl)
        y0 = convert(Vector{T},y0)
        return ActivityBubbleTemperature{T}(vol0,T0,y0,nonvolatiles,itmax_ss,rtol_ss,gas_fug,poynting)
    elseif  !isnothing(T0) && !isnothing(y0)
        T0,_ = promote(T0,first(y0))
        T = eltype(T0)
        y0 = convert(Vector{T},y0)
        return ActivityBubbleTemperature{T}(vol0,T0,y0,nonvolatiles,itmax_ss,rtol_ss,gas_fug,poynting)
    else
        throw(error("invalid specification for bubble temperature"))
    end
end

export ActivityBubblePressure,ActivityBubbleTemperature