#this model only holds a named tuple with all models.

"""
    CompositeModel(components;
    gas = BasicIdeal,
    liquid = RackettLiquid,
    saturation = LeeKeslerSat,
    gas_userlocations = String[]
    liquid_userlocations = String[]
    saturation_userlocations = String[]

Composite Model. it is not consistent, but it can hold different correlations that
are faster than a volume or saturation pressure iteration.

"""
struct CompositeModel{NT} <: EoSModel
    components::Vector{String}
    models::NT
end

function volume(model::CompositeModel,p,T,z=SA[1.0];phase=:unknown,threaded=false)
    if is_liquid(phase)
        return volume(model.models.liquid,p,T,z;phase,threaded)
    elseif is_vapour(phase)
        return volume(model.models.gas,p,T,z;phase,threaded)
    else
        if length(model) == 1
            psat,vl,vv = saturation_pressure(model,T)
            if !isnan(psat)
                if p > psat
                    return vl
                else
                    return vv
                end
            else
                tc,pc,vc = crit_pure(model)
                if T > tc #supercritical conditions. ideally, we could go along the critical isochore, but we dont have that.
                    if p > pc # supercritical fluid
                        return volume(model.models.liquid,p,T,z;phase,threaded)
                    else #gas phase
                        return volume(model.models.gas,p,T,z;phase,threaded)
                    end
                else #something failed on saturation_pressure, not related to passing the critical point
                    @error "an error ocurred while determining saturation line division."
                    return nan
                end

            end
            @error "A phase needs to be specified on multicomponent composite models."
            return nan
        end
    end
end

function saturation_pressure(model::CompositeModel,T,v0 = nothing)
    nan = zero(T)/zero(T)
    psat,_,_ = saturation_pressure(model.saturation,T)
    if !isnan(psat)
        vl = volume(model.liquid,psat,T,phase=:l)
        vv = volume(model.gas,psat,T,phase=:v)
        return psat,vl,vv
    #if psat fails, there are two options:
    #1- over critical point -> nan nan nan
    #2- saturation failed -> nan nan nan
    else 
        return nan,nan,nan
    end
end

function crit_pure(model::CompositeModel)
    return crit_pure(model.saturation)
end

include("SaturationModel/SaturationModel.jl")



