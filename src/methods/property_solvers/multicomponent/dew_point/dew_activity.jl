"""
    ActivityDewPressure(kwargs...)

Function to compute [`dew_pressure`](@ref) using Activity Coefficients.
On activity coefficient models it solves the problem via succesive substitucion.
On helmholtz-based models, it uses the Chapman approximation for activity coefficients.

Inputs:
- `gas_fug = true`: if the solver uses gas fugacity coefficients. on `ActivityModel` is set by default to `false`
- `poynting = true`: if the solver use the poynting correction on the liquid fugacity coefficients. on `ActivityModel` is set by default to `false`
- `x0 = nothing`: optional, initial guess for the liquid phase composition
- `p0 = nothing`: optional, initial guess for the dew pressure `[Pa]`
- `vol0 = nothing`: optional, initial guesses for the liquid and vapor phase volumes `[m³]`
- `atol = 1e-8`: optional, absolute tolerance of the non linear system of equations
- `rtol = 1e-12`: optional, relative tolerance of the non linear system of equations
- `itmax_ss = 40`: optional, maximum number of sucesive substitution iterations
"""
struct ActivityDewPressure{T} <: DewPointMethod
    vol0::Union{Nothing,Tuple{T,T}}
    p0::Union{Nothing,T}
    x0::Union{Nothing,Vector{T}}
    nonvolatiles::Union{Nothing,Vector{String}}
    itmax_ss::Int64
    rtol_ss::Float64
    gas_fug::Bool
    poynting::Bool
end

function ActivityDewPressure(;vol0 = nothing,
                                p0 = nothing,
                                x0 = nothing,
                                nonvolatiles = nothing,
                                itmax_ss = 40,
                                rtol_ss = 1e-8,
                                gas_fug = true,
                                poynting = true)

    if p0 == x0 == vol0 == nothing
        return ActivityDewPressure{Nothing}(vol0,p0,x0,nonvolatiles,itmax_ss,rtol_ss,gas_fug,poynting)
    elseif (p0 == x0 == nothing) && !isnothing(vol0)
        vl,vv = promote(vol0[1],vol0[2])
        return ActivityDewPressure{typeof(vl)}(vol0,p0,x0,nonvolatiles,itmax_ss,rtol_ss,gas_fug,poynting)
    elseif (vol0 == x0 == nothing) && !isnothing(p0)
        p0 = float(p0)
        return ActivityDewPressure{typeof(p0)}(vol0,p0,x0,nonvolatiles,itmax_ss,rtol_ss,gas_fug,poynting)
    elseif (p0 == vol0 == nothing) && !isnothing(x0)
        T = eltype(x0)
        return ActivityDewPressure{T}(vol0,p0,x0,nonvolatiles,itmax_ss,rtol_ss,gas_fug,poynting)
    elseif !isnothing(vol0) && !isnothing(p0) && !isnothing(x0)
        vl,vv,p0,_ = promote(vol0[1],vol0[2],p0,first(x0))
        T = eltype(vl)
        x0 = convert(Vector{T},x0)
        return ActivityDewPressure{T}(vol0,p0,x0,nonvolatiles,itmax_ss,rtol_ss,gas_fug,poynting)
    elseif !isnothing(vol0) && !isnothing(x0)
        vl,vv,_ = promote(vol0[1],vol0[2],first(x0))
        T = eltype(vl)
        x0 = convert(Vector{T},x0)
        return ActivityDewPressure{T}(vol0,p0,x0,nonvolatiles,itmax_ss,rtol_ss,gas_fug,poynting)
    elseif  !isnothing(p0) && !isnothing(x0)
        p0,_ = promote(p0,first(x0))
        T = eltype(p0)
        x0 = convert(Vector{T},x0)
        return ActivityDewPressure{T}(vol0,p0,x0,nonvolatiles,itmax_ss,rtol_ss,gas_fug,poynting)
    else
        throw(error("invalid specification for dew pressure"))
    end
end


function dew_pressure_impl(model,T,y,method::ActivityDewPressure)
    R̄ = Rgas(model)
    pure = split_pure_model(model)
    sat = saturation_pressure.(pure,T)
    vl_pure = getindex.(sat,2)
    p_pure = first.(sat)
    vv_pure = last.(sat)
    if isnothing(method.vol0)
        vv = dot(vv_pure,y)
        vl = zero(vv)
    else
        vl,vv = method.vol0
    end
   
    if isnan(vv)
        return vv,vv,vv,y
    end

    if isnothing(method.p0)
        pmix = pressure(model,vv,T,y)
    else
        pmix = method.p0
    end

    μmix = zeros(typeof(pmix),length(pure))
    ϕ = copy(μmix)
    logϕ = copy(μmix)
    x = copy(μmix)
    ϕpure = copy(μmix)
    ϕpure .= 1
    x .= y .* pmix ./ p_pure
    x ./= sum(x)
    
    if !isnothing(method.x0)
        x .= method.x0
    end

    if iszero(vl)
        vl = dot(vl_pure,x)
    end

    μmix = VT_chemical_potential_res!(μmix,model,vl,T,x)
    RT = (R̄*T)
    ϕ .= 1
    μpure = only.(VT_chemical_potential_res.(pure,vl_pure,T))
    if method.gas_fug
        if iszero(vv)
            vv = dot(last.(sat),y)
        end
        μv = VT_chemical_potential_res!(ϕ,model,vv,T,y)
        ϕ .= exp.(μv ./ RT .- log.(pmix .* vv ./ RT))
        ϕpure .= exp.(μpure ./ RT .- log.(p_pure .* vl_pure ./ RT))
    end
    pold = zero(pmix)
    γ = zeros(typeof(pmix),length(pure))
    #pure part
    
    
    use_𝒫 = method.poynting
    if use_𝒫
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
            if use_𝒫
                ln𝒫 = vpureᵢ*expm1(κ[i]*(pmix-pᵢ))/(κ[i]*RT) #see end of file
                𝒫 = exp(ln𝒫)
            else
                𝒫 = one(pᵢ)
            end
            
            #y[i]*ϕ[i]*P = x[i]*γ[i]*pᵢ*ϕ̂ᵢ*𝒫
            x[i] = y[i]*ϕ[i]/(γ[i]*pᵢ*ϕ̂ᵢ*𝒫)
            #y[i] = x[i]*γ[i]*pᵢ*𝒫*ϕ̂ᵢ/ϕ[i]
        end
        pold = pmix
        pmix = 1/sum(x)
        
        x ./= sum(x)
        vl = volume(model,pmix,T,x,vol0 = vl)
        if method.gas_fug
            logϕ, vv = lnϕ!(logϕ,model,pmix,T,y,phase = :vapor, vol0 = vv)
            ϕ .= exp.(logϕ)
        else
            vv = volume(model,pmix,T,y,phase =:vapor,vol0 = vv)
        end
        err = abs(pold-pmix)/pmix
        μmix = VT_chemical_potential_res!(μmix,model,vl,T,x)
        if err < method.rtol_ss
            break
        end
        if !isfinite(err)
            break
        end
    end
    return pmix,vl,vv,x
end

"""
    ActivityDewTemperature(kwargs...)

Function to compute [`dew_temperature`](@ref) using Activity Coefficients.
On activity coefficient models it solves the problem via succesive substitucion.
On helmholtz-based models, it uses the Chapman approximation for activity coefficients.

Inputs:
- `gas_fug = true`: if the solver uses gas fugacity coefficients. on `ActivityModel` is set by default to `false`
- `poynting = true`: if the solver use the poynting correction on the liquid fugacity coefficients. on `ActivityModel` is set by default to `false`
- `x0 = nothing`: optional, initial guess for the liquid phase composition
- `T0 = nothing`: optional, initial guess for the dew temperature `[K]`
- `vol0 = nothing`: optional, initial guesses for the liquid and vapor phase volumes `[m³]`
- `atol = 1e-8`: optional, absolute tolerance of the non linear system of equations
- `rtol = 1e-12`: optional, relative tolerance of the non linear system of equations
- `itmax_ss = 40`: optional, maximum number of sucesive substitution iterations
"""
struct ActivityDewTemperature{T} <: DewPointMethod
    vol0::Union{Nothing,Tuple{T,T}}
    T0::Union{Nothing,T}
    x0::Union{Nothing,Vector{T}}
    nonvolatiles::Union{Nothing,Vector{String}}
    itmax_ss::Int64
    rtol_ss::Float64
    gas_fug::Bool
    poynting::Bool
end

function ActivityDewTemperature(;vol0 = nothing,
                                T0 = nothing,
                                x0 = nothing,
                                nonvolatiles = nothing,
                                itmax_ss = 40,
                                rtol_ss = 1e-8,
                                gas_fug = true,
                                poynting = true)

    if T0 == x0 == vol0 == nothing
        return ActivityDewTemperature{Nothing}(vol0,T0,x0,nonvolatiles,itmax_ss,rtol_ss,gas_fug,poynting)
    elseif (T0 == x0 == nothing) && !isnothing(vol0)
        vl,vv = promote(vol0[1],vol0[2])
        return ActivityDewTemperature{typeof(vl)}(vol0,T0,x0,nonvolatiles,itmax_ss,rtol_ss,gas_fug,poynting)
    elseif (vol0 == x0 == nothing) && !isnothing(T0)
        T0 = float(T0)
        return ActivityDewTemperature{typeof(T0)}(vol0,T0,x0,nonvolatiles,itmax_ss,rtol_ss,gas_fug,poynting)
    elseif (T0 == vol0 == nothing) && !isnothing(T0)
        T = eltype(x0)
        return ActivityDewTemperature{T}(vol0,T0,x0,nonvolatiles,itmax_ss,rtol_ss,gas_fug,poynting)
    elseif !isnothing(vol0) && !isnothing(T0) && !isnothing(x0)
        vl,vv,T0,_ = promote(vol0[1],vol0[2],T0,first(x0))
        T = eltype(vl)
        x0 = convert(Vector{T},x0)
        return ActivityDewTemperature{T}(vol0,T0,x0,nonvolatiles,itmax_ss,rtol_ss,gas_fug,poynting)
    elseif !isnothing(vol0) && !isnothing(x0)
        vl,vv,_ = promote(vol0[1],vol0[2],first(x0))
        T = eltype(vl)
        x0 = convert(Vector{T},x0)
        return ActivityDewTemperature{T}(vol0,T0,x0,nonvolatiles,itmax_ss,rtol_ss,gas_fug,poynting)
    elseif  !isnothing(T0) && !isnothing(x0)
        T0,_ = promote(T0,first(x0))
        T = eltype(T0)
        x0 = convert(Vector{T},x0)
        return ActivityDewTemperature{T}(vol0,T0,x0,nonvolatiles,itmax_ss,rtol_ss,gas_fug,poynting)
    else
        throw(error("invalid specification for dew temperature"))
    end
end

export ActivityDewPressure,ActivityDewTemperature