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
        throw(error("invalid specification for bubble pressure"))
    end
end


function dew_pressure_impl(model,T,y,method::ActivityDewPressure)
    
    pure = split_model(model)
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
        return vv,vv,vv,x
    end

    if isnothing(method.p0)
        pmix = pressure(model,vv,T,y)
    else
        pmix = method.p0
    end

    μmix = zeros(typeof(pmix),length(pure))
    ϕ = copy(μmix)
    x = copy(μmix)
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
    if method.gas_fug
        if iszero(vv)
            vv = dot(last.(sat),y)
        end
        μv = VT_chemical_potential_res!(ϕ,model,vv,T,y)
        ϕ .= exp.(μv ./ RT .- log.(pmix .* vv ./ RT))
    end
    pold = zero(pmix)
    γ = zeros(typeof(pmix),length(pure))
    #pure part
    μpure = only.(VT_chemical_potential_res.(pure,vl_pure,T))
    ϕpure = exp.(μpure ./ RT .- log.(p_pure .* vl_pure ./ RT))
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
            ϕ̂ᵢ =  ϕpure[i]
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
            logϕ, vv = lnϕ(model,pmix,T,y,phase = :vapor, vol0 = vv)
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