function dew_pressure(model::ActivityModel,T,x)
    dew_pressure(model,T,x,ActivityDewPressure(gas_fug = false, poynting = false))
end

function dew_pressure(model::ActivityModel, T, y, method::DewPointMethod)
    if !(method isa ActivityDewPressure)
        throw(error("$method not supported by Activity models"))
    end
    _T = typeof(T)
    _y = typeof(y)
    Base.invoke(dew_pressure,Tuple{EoSModel,_T,_y,DewPointMethod},model,T,y,method)
end

function dew_pressure_impl(model::ActivityModel,T,y,method::ActivityDewPressure)
    sat = saturation_pressure.(model.puremodel,T)
    pmodel = model.puremodel.model
    pure = model.puremodel.pure
    p_pure = first.(sat)
    vl_pure = getindex.(sat,2)
    vv_pure = last.(sat)
    vv = dot(vv_pure,y)
    x = y ./ p_pure #raoult
    p = 1/sum(x)
    x .*= p
    if isnan(p)
        return (p,p,p,x)
    end

    vl = volume(pmodel,p,T,x,phase = :l)
    γ = activity_coefficient(model,1e-4,T,x)
    ϕ = copy(x)
    ϕ .= 1.0
    RT = R̄*T
    if method.gas_fug
        logϕ, vv = lnϕ(pmodel,p,T,y,phase = :vapor, vol0 = vv)
        ϕ .= exp.(logϕ)
    else
        vv = volume(pmodel,p,T,y,phase = :vapor, vol0 = vv)
    end
    #fugacity of pure component at saturation conditions
    if method.gas_fug
        μpure = only.(VT_chemical_potential_res.(pure,vv_pure,T))
        ϕpure = exp.(μpure ./ RT .- log.(p_pure .* vv_pure ./ RT))
    else
        ϕpure = copy(ϕ)
        ϕpure .= 1.0
    end
    pold = zero(p)
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
            x[i] = y[i]*ϕ[i]/(γ[i]*pᵢ*ϕ̂ᵢ*𝒫)
        end
        pold = p
        p = 1/sum(x)
        x .*= p
        if method.gas_fug
            logϕ, vv = lnϕ(pmodel,p,T,y,phase = :vapor, vol0 = vv)
            ϕ .= exp.(logϕ)
        else
            vv = volume(pmodel,p,T,y,phase = :vapor, vol0 = vv)
        end
        err = abs(pold-p)/p
        γ .= activity_coefficient(model,p,T,x)
        !isfinite(err) && break
        err < method.rtol_ss && break
    end
    vl = volume(pmodel,p,T,x,phase = :liquid,vol0 = vl)
    return (p,vl,vv,x)
end