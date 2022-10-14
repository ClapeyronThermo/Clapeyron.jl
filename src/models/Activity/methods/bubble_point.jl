function bubble_pressure(model::ActivityModel,T,x)
    bubble_pressure(model,T,x,ActivityBubblePressure(gas_fug = false, poynting = false))
end

function bubble_pressure(model::ActivityModel, T, x, method::BubblePointMethod)
    if !(method isa ActivityBubblePressure)
        throw(error("$method not supported by Activity models"))
    end
    _T = typeof(T)
    _x = typeof(x)
    Base.invoke(bubble_pressure,Tuple{EoSModel,_T,_x,BubblePointMethod},model,T,x,method)
end


function bubble_pressure_impl(model::ActivityModel,T,x,method::ActivityBubblePressure)
    sat = saturation_pressure.(model.puremodel,T)
    pmodel = model.puremodel.model
    pure = model.puremodel.pure
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
    ϕ .= 1.0
    RT = R̄*T
    if method.gas_fug
        logϕ, vv = lnϕ(__gas_model(pmodel),p,T,y,phase = :vapor, vol0 = vv)
        ϕ .= exp.(logϕ)
    else
        vv = volume(pmodel,p,T,y,phase = :vapor, vol0 = vv)
    end
    #fugacity of pure component at saturation conditions
    if method.gas_fug
        μpure = only.(VT_chemical_potential_res.(__gas_model.(pure),vv_pure,T))
        ϕpure = exp.(μpure ./ RT .- log.(p_pure .* vv_pure ./ RT))
    else
        ϕpure = copy(ϕ)
        ϕpure .= 1.0
    end
    for k in 1:method.itmax_ss
        for i in eachindex(γ)
            pᵢ = p_pure[i]
            vpureᵢ = vl_pure[i]
            ϕ̂ᵢ =  ϕpure[i]
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
            logϕ, vv = lnϕ(__gas_model(pmodel),p,T,y,phase = :vapor, vol0 = vv)
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

function bubble_temperature(model::ActivityModel,T,x)
    bubble_temperature(model,T,x,ActivityBubbleTemperature(gas_fug = false, poynting = false))
end

function bubble_temperature(model::ActivityModel, T, x, method::BubblePointMethod)
    if !(method isa ActivityBubbleTemperature)
        throw(error("$method not supported by Activity models"))
    end
    _T = typeof(T)
    _x = typeof(x)
    Base.invoke(bubble_temperature,Tuple{EoSModel,_T,_x,BubblePointMethod},model,T,x,method)
end

function bubble_temperature_impl(model::ActivityModel,p,x,method::ActivityBubbleTemperature)
    f(z) = Obj_bubble_temperature(model,z,p,x)
        pure = model.puremodel
        sat = saturation_temperature.(pure,p)
        Ti   = first.(sat)
        T0 = dot(Ti,x)
        T = Roots.find_zero(f,T0)
    p,vl,vv,y = bubble_pressure(model,T,x,ActivityBubblePressure(gas_fug = method.gas_fug,poynting = method.poynting))
    return (T,vl,vv,y)
end

function Obj_bubble_temperature(model::ActivityModel,T,p,x)
    sat = saturation_pressure.(model.puremodel,T)
    p_sat = [tup[1] for tup in sat]
    γ     = activity_coefficient(model,p,T,x)
    y     = x.*γ.*p_sat ./ p
    return sum(y)-1
end
