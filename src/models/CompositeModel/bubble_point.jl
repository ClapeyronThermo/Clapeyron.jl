function init_preferred_method(method::typeof(bubble_pressure),model::RestrictedEquilibriaModel,kwargs)
    gas_fug = get(kwargs,:gas_fug,false)
    poynting = get(kwargs,:poynting,false)
    return ActivityBubblePressure(;gas_fug,poynting,kwargs...)
end

function bubble_pressure_impl(model::RestrictedEquilibriaModel,T,x,method::ActivityBubblePressure) 
    if model isa GammaPhi
        pmodel = model.fluid.model
        pure = model.fluid.pure
    else
        pmodel = model
        pure = split_model(pmodel,1:length(model))
    end

    sat = saturation_pressure.(pure,T)
    p_pure = first.(sat)
    vl_pure = getindex.(sat,2)
    vv_pure = last.(sat)
    γ = activity_coefficient(model,1e-4,T,x)
    y = x.*γ.*p_pure
    p = sum(y)
    y ./= p
    if isnan(p)
        return (p,p,p,y)
    end
    vv = volume(pmodel,p,T,y,phase = :v)
    if !method.gas_fug && !method.poynting
        vl = volume(pmodel,p,T,x,phase = :l)
        return (p,vl,vv,y)
    end

    ϕ = copy(y)
    logϕ = copy(y)
    ϕ .= 1.0
    RT = R̄*T
    if method.gas_fug
        logϕ, vv = lnϕ!(logϕ,gas_model(pmodel),p,T,y,phase = :vapor, vol0 = vv)
        ϕ .= exp.(logϕ)
    else
        vv = volume(pmodel,p,T,y,phase = :vapor, vol0 = vv)
    end
    #fugacity of pure component at saturation conditions
    if method.gas_fug
        μpure = only.(VT_chemical_potential_res.(gas_model.(pure),vv_pure,T))
        ϕpure = exp.(μpure ./ RT .- log.(p_pure .* vv_pure ./ RT))
    else
        ϕpure = copy(ϕ)
        ϕpure .= 1.0
    end
    for k in 1:method.itmax_ss
        for i in eachindex(γ)
            pᵢ = p_pure[i]
            vpureᵢ = vl_pure[i]
            ϕ̂ᵢ = ϕpure[i]
            if method.poynting && method.gas_fug
                ln𝒫 = vpureᵢ*(p - pᵢ)/RT
                𝒫 = exp(ln𝒫)
            else
                𝒫 = one(pᵢ)
            end
            y[i] = x[i]*γ[i]*pᵢ*𝒫*ϕ̂ᵢ/ϕ[i] #really yᵢ*P, we normalize later
        end
        pold = p
        p = sum(y)
        y ./= p
        if method.gas_fug
            logϕ, vv = lnϕ!(logϕ,gas_model(pmodel),p,T,y,phase = :vapor, vol0 = vv)
            ϕ .= exp.(logϕ)
        else
            vv = volume(pmodel,p,T,y,phase = :vapor, vol0 = vv)
        end
        err = abs(pold-p)/p
        !isfinite(err) && break
        err < method.rtol_ss && break
    end
    vl = volume(pmodel,p,T,x,phase = :l)
    return (p,vl,vv,y)
end

function init_preferred_method(method::typeof(bubble_temperature),model::RestrictedEquilibriaModel,kwargs)
    gas_fug = get(kwargs,:gas_fug,false)
    poynting = get(kwargs,:poynting,false)
    return ActivityBubbleTemperature(;gas_fug,poynting,kwargs...)
end

function bubble_temperature_impl(model::RestrictedEquilibriaModel,p,x,method::ActivityBubbleTemperature)
    if model isa GammaPhi
        pmodel = model.fluid.model
        pure_all = model.fluid.pure
        if model.fluid.model isa FluidCorrelation
            pure = map(x -> x.saturation,pure_all)
        else
            pure = pure_all
        end
    else
        pmodel = model
        pure = split_model(pmodel,1:length(model))
    end

    f(z) = Obj_bubble_temperature(model,z,p,x,pure)
    sat = saturation_temperature.(pure,p)
    Ti   = first.(sat)
    T0 = dot(Ti,x)
    T = Roots.find_zero(f,T0)
    p,vl,vv,y = bubble_pressure(model,T,x,ActivityBubblePressure(gas_fug = method.gas_fug,poynting = method.poynting))
    return (T,vl,vv,y)
end

function Obj_bubble_temperature(model::GammaPhi,T,p,x,pure)
    sat = saturation_pressure.(pure,T)
    p_sat = [tup[1] for tup in sat]
    γ     = activity_coefficient(model,p,T,x)
    y     = x.*γ.*p_sat ./ p
    return sum(y)-1
end
