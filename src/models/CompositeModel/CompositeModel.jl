#this model only holds a named tuple with all models.
include("SaturationModel/SaturationModel.jl")
include("LiquidVolumeModel/LiquidVolumeModel.jl")
"""
    CompositeModel(components;
    gas = BasicIdeal,
    liquid = RackettLiquid,
    saturation = LeeKeslerSat,
    gas_userlocations = String[],
    liquid_userlocations = String[],
    saturation_userlocations = String[]

Composite Model. it is not consistent, but it can hold different correlations that
are faster than a volume or saturation pressure iteration.

"""
struct CompositeModel{𝕍,𝕃,𝕊,𝕃𝕍,𝕃𝕊} <: EoSModel
    components::Vector{String}
    gas::𝕍
    liquid::𝕃
    solid::𝕊
    saturation::𝕃𝕍
    melting::𝕃𝕊
end

Base.length(cmodel::CompositeModel) = length(cmodel.components)

function CompositeModel(components;
    liquid = RackettLiquid,
    gas = BasicIdeal,
    solid = nothing,
    saturation = LeeKeslerSat,
    melting = nothing,
    gas_userlocations = String[],
    liquid_userlocations = String[],
    solid_userlocations = String[],
    saturation_userlocations = String[],
    melting_userlocations = String[],
    verbose = false)

    init_gas = init_model(gas,components,gas_userlocations,verbose)
    init_liquid = init_model(liquid,components,liquid_userlocations,verbose)
    init_solid = init_model(solid,components,solid_userlocations,verbose)
    init_sat = init_model(saturation,components,saturation_userlocations,verbose)
    init_melt = init_model(melting,components,melting_userlocations,verbose)
    return CompositeModel(components,init_gas,init_liquid,init_solid,init_sat,init_melt)
end

function Base.show(io::IO,mime::MIME"text/plain",model::CompositeModel)
    print(io,"Composite Model:")
    model.gas !== nothing && print(io,'\n'," Gas Model: ",model.gas)
    model.liquid !== nothing && print(io,'\n'," Liquid Model: ",model.liquid)
    model.solid !== nothing && println(io,'\n'," Solid Model: ",model.solid)
    model.saturation !== nothing && print(io,'\n'," Saturation Model: ",model.saturation)
    model.melting !== nothing && print(io,'\n'," Melting Model: ",model.melting)
end

function Base.show(io::IO,model::CompositeModel)
    eosshow(io,model)
end

__gas_model(model::CompositeModel) = model.gas

function volume_impl(model::CompositeModel,p,T,z,phase=:unknown,threaded=false,vol = vol0)
    _0 = zero(p+T+first(z))
    nan = _0/_0
    if is_liquid(phase)
        return volume(model.liquid,p,T,z;phase,threaded)
    elseif is_vapour(phase)
        return volume(model.gas,p,T,z;phase,threaded)
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
                        return volume(model.liquid,p,T,z;phase,threaded)
                    else #gas phase
                        return volume(model.gas,p,T,z;phase,threaded)
                    end
                else #something failed on saturation_pressure, not related to passing the critical point
                    @error "an error ocurred while determining saturation line division."
                    return nan
                end

            end
            
            return nan
        else
            @error "A phase needs to be specified on multicomponent composite models."
            return nan
        end
    end
end

function saturation_pressure(model::CompositeModel,T::Real)
    if model.saturation isa SaturationModel
        method = SaturationCorrelation()
    else
        method = ChemPotVSaturation()
    end
    return saturation_pressure(model,T,method)
end

function saturation_pressure(model::CompositeModel,T,method::SaturationMethod)
    nan = zero(T)/zero(T)
    psat,_,_ = saturation_pressure(model.saturation,T,method)
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
    return crit_pure(model.models.saturation)
end

function saturation_temperature(model::CompositeModel,p,method::SaturationMethod)
    nan = zero(p)/zero(p)
    Tsat,_,_ = saturation_temperature(model.saturation,p,method)
    if !isnan(Tsat)
        vl = volume(model.liquid,p,Tsat,phase=:l)
        vv = volume(model.gas,p,Tsat,phase=:v)
        return Tsat,vl,vv
    #if psat fails, there are two options:
    #1- over critical point -> nan nan nan
    #2- saturation failed -> nan nan nan
    else 
        return nan,nan,nan
    end
end

export CompositeModel
