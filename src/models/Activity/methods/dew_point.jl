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
            logϕ, vv = lnϕ(__gas_model(pmodel),p,T,y,phase = :vapor, vol0 = vv)
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

function dew_temperature(model::ActivityModel,T,x)
    dew_temperature(model,T,x,ActivityDewTemperature(gas_fug = false, poynting = false))
end

function dew_temperature(model::ActivityModel, T, x, method::DewPointMethod)
    if !(method isa ActivityDewTemperature)
        throw(error("$method not supported by Activity models"))
    end
    _T = typeof(T)
    _x = typeof(x)
    Base.invoke(dew_temperature,Tuple{EoSModel,_T,_x,DewPointMethod},model,T,x,method)
end

function dew_temperature_impl(model::ActivityModel,p,y,method::ActivityDewTemperature)
    f(z) = Obj_bubble_temperature(model,z,p,y)
    pure = model.puremodel
    sat = saturation_temperature.(pure,p)
    Ti   = first.(sat)
    T0 = dot(Ti,y)
    sat0 = saturation_pressure.(pure,T0)
    pi0 = first.(sat0)
    x0 = y ./ pi0
    x0 ./= sum(x0)
    x0[end] = T0
    f0(F,w) =  Obj_dew_temperature(F,model,p,y,w[1:end-1],w[end])
    sol = Solvers.nlsolve(f0,x0,LineSearch(Newton()))
    wsol = Solvers.x_sol(sol)
    T = wsol[end]
   
    x = collect(FractionVector(wsol[1:end-1]))
    vl = volume(pure.model,p,T,x,phase = :l)
    vv = volume(pure.model,p,T,y,phase = :v)
    return (T,vl,vv,x)
    #p,vl,vv,y = bubble_pressure(model,T,x,ActivityBubblePressure(gas_fug = method.gas_fug,poynting = method.poynting))
    #return (T,vl,vv,y)
end

function Obj_dew_temperature(F,model::ActivityModel,p,y,_x,T)
    x = FractionVector(_x)
    γ  = activity_coefficient(model,p,T,x)
    sat = saturation_pressure.(model.puremodel,T) #TODO: AD rule for saturation pressure
    for i in eachindex(F)
        pᵢ = sat[i][1]
        F[i] = x[i] - y[i]*p/(γ[i]*pᵢ)
    end
    return F
end
