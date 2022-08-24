
function bubble_pressure(model::ActivityModel,T,x)
    bubble_pressure(model,T,x,ActivityBubblePressure(gas_fug = false, poynting = false))
end

function bubble_pressure_impl(model::ActivityModel,T,x,method::BubblePointMethod)
    throw(error("$method not supported by Activity models"))
end

function bubble_pressure_impl(model::ActivityModel,T,x,method::ActivityBubblePressure)
    sat = saturation_pressure.(model.puremodel,T)
    pmodel = model.puremodel.model
    pure = model.puremodel.pure
    p_pure = first.(sat)
    vl_pure = getindex.(sat,2)
    vv_pure = last.(sat)
    γ     = activity_coefficient(model,1e-4,T,x)
    p     = sum(x.*γ.*p_pure)
    y     = x.*γ.*p_pure ./ p
    vl = volume(pmodel,p,T,x,phase = :l)
    vv = volume(pmodel,p,T,y,phase = :v)
    if !method.gas_fug && !method.poynting
        return (p,vl,vv,y)
    end
    
    ϕ = copy(y)
    ϕ .= 1.0
    RT = R̄*T
    if method.gas_fug
        logϕ, vv = lnϕ(pmodel,p,T,y,phase = :vapor, vol0 = vv)
        ϕ .= exp.(logϕ)
    else
        vv = volume(pmodel,p,T,y,phase = :vapor, vol0 = vv)
    end
    μpure = only.(VT_chemical_potential_res.(pure,vl_pure,T))
    ϕpure = exp.(μpure ./ RT .- log.(p_pure .* vl_pure ./ RT))
    for k in 1:method.itmax_ss
        for i in eachindex(γ)
            pᵢ = p_pure[i]
            vpureᵢ = vl_pure[i]
            ϕ̂ᵢ =  ϕpure[i]
            if method.poynting
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
        vl = volume(pmodel,p,T,x,phase = :liquid,vol0 = vl)
        if method.gas_fug
            logϕ, vv = lnϕ(pmodel,p,T,y,phase = :vapor, vol0 = vv)
            ϕ .= exp.(logϕ)
        else
            vv = volume(pmodel,p,T,y,phase = :vapor, vol0 = vv)
        end
        err = abs(pold-p)/p
        if err < method.rtol_ss
            break
        end
    end
    return (p,vl,vv,y)
end

function bubble_temperature(model::ActivityModel,p,x)
    f(z) = Obj_bubble_temperature(model,z,p,x)
    if T0===nothing
        pure = model.puremodel
        sat = saturation_temperature.(pure,p)
        Ti   = zero(x)
        for i ∈ 1:length(x)
            if isnan(sat[i][1])
                Tc,pc,vc = crit_pure(pure[i])
                g(x) = p-pressure(pure[i],vc,x,[1.])
                Ti[i] = Roots.find_zero(g,(Tc))
            else
                Ti[i] = sat[i][1]
            end
        end
        T = Roots.find_zero(f,(minimum(Ti)*0.9,maximum(Ti)*1.1))
    else
        T = Roots.find_zero(f,T0)
    end
    p,vl,vv,y = bubble_pressure(model,T,x)
    return (T,vl,vv,y)
end

function Obj_bubble_temperature(model::ActivityModel,T,p,x)
    sat = saturation_pressure.(model.puremodel,T)
    p_sat = [tup[1] for tup in sat]
    γ     = activity_coefficient(model,1e-4,T,x)
    y     = x.*γ.*p_sat ./ p
    return sum(y)-1
end
