struct ActivityBubblePressure{T,M} <: BubblePointMethod
    y0::Union{Vector{T},Nothing}
    switch_method::M
end

function bubble_pressure_chapman(model,T,x;v0 = nothing)
    pure = split_model(model)
    sat = saturation_pressure.(pure,T)
    vl_pure = getindex.(sat,2)
    p_pure = first.(sat)
    vv_pure = last.(sat)
    if isnothing(v0)
        vmix = dot(vl_pure,x)
    else
        vmix = v0
    end
    vgas = zero(vmix)
    if isnan(vmix)
        return vmix,vmix,vmix,x
    end
    μmix = VT_chemical_potential_res(model,vmix,T,x)
    ϕ = copy(μmix)
    ϕ .= 1
    pmix = pressure(model,vmix,T,x)
    pold = zero(pmix)
    γ = zeros(typeof(pmix),length(pure))
    y = zeros(length(pure))
    RT = (R̄*T)
    μpure = only.(VT_chemical_potential_res.(pure,vl_pure,T))
    ϕpure = exp.(μpure ./ RT .- log.(p_pure .* vl_pure ./ RT))
    κ = VT_isothermal_compressibility.(pure,vl_pure,T)
    for k in 1:20
        for i in eachindex(γ)
            pᵢ = p_pure[i]
            vpureᵢ = vl_pure[i]
            μᵢ = μpure[i]
            ϕ̂ᵢ =  ϕpure[i]
            γ[i] = exp(log(vpureᵢ/vmix) + (μmix[i] - μᵢ)/RT -  vpureᵢ*(pmix -pᵢ)/RT)
            #𝒫 = exp(vpureᵢ*(pmix-pᵢ)/RT)
            #ln𝒫1 = vpureᵢ*(pmix-pᵢ)/RT# - dvᵢ∂P[i]*(0.5*(pmix^2 - pᵢ^2) - pᵢ*(pmix-pᵢ))/RT
            ln𝒫 = vpureᵢ*expm1(κ[i]*(pmix-pᵢ))/(κ[i]*RT)
            𝒫 = exp(ln𝒫)
            #𝒫 = ifelse(𝒫_chapman > 2,𝒫_simple,𝒫_chapman)
            y[i] = x[i]*γ[i]*pᵢ*𝒫*ϕ̂ᵢ/ϕ[i]
        end
        pold = pmix
        pmix = sum(y)
        y ./= pmix
        if iszero(vgas)
            vgas = dot(y,vv_pure)
        end
        logϕ, vgas = lnϕ(model,pmix,T,y,phase = :vapor, vol0 = vgas)
        vmix = volume(model,pmix,T,x,vol0 = vmix)
        ϕ .= exp.(logϕ)
        err = abs(pold-pmix)
        μmix = VT_chemical_potential_res!(μmix,model,vmix,T,x)
        if err < 2e-8
            break
        end
        @show err
    end
    
    return pold,vmix,vgas,y
end
## Poynting factor, Clapeyron edition
## we calculate ∫vi/RT dP from Ppure to Pmix
## v(p) = v0*exp(Δ(p)) = v0*exp((p-p0)*κ) #v0, p0 = vi, pi
## v(p) = vi*exp(-pi*κ) * exp(κ*p)
## ∫v(p)/RT dp = vi*exp(-pi*κ)/RT*∫exp(κ*p) dp = vi*exp(-pi*κ)*exp(κ*p)/κ 
## vi*exp(κ*(p-pi))/κRT  |from pi to pmix

