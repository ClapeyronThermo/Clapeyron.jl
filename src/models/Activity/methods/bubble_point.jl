function bubble_pressure(model::ActivityModel,T,x)
    sat = saturation_pressure.(model.puremodel,T)
    p_sat = [tup[1] for tup in sat]
    γ     = activity_coefficient(model,1e-4,T,x)
    p     = sum(x.*γ.*p_sat)
    y     = x.*γ.*p_sat ./ p
    vl = volume(model.puremodel.model,p,T,x,phase = :l)
    vv = volume(model.puremodel.model,p,T,y,phase = :v)
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
