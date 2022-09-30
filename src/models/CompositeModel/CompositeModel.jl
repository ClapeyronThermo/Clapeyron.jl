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
struct CompositeModel{NT} <: EoSModel
    components::Vector{String}
    models::NT
end

Base.length(cmodel::CompositeModel) = length(cmodel.components)

function CompositeModel(components;
    gas = BasicIdeal,
    liquid = RackettLiquid,
    saturation = LeeKeslerSat,
    gas_userlocations = String[],
    liquid_userlocations = String[],
    saturation_userlocations = String[],
    verbose = false)

    init_gas = init_model(gas,components,gas_userlocations,verbose)
    init_liquid = init_model(liquid,components,liquid_userlocations,verbose)
    init_sat = init_model(saturation,components,saturation_userlocations,verbose)
    models = (gas = init_gas,liquid = init_liquid,saturation = init_sat)
    return CompositeModel(components,models)
end

function Base.show(io::IO,mime::MIME"text/plain",model::CompositeModel)
    println(io,"Composite Model:")
    println(io," Liquid Model: ",model.models.liquid)
    println(io," Gas Model: ",model.models.gas)
    print(io," Saturation Model: ",model.models.saturation)
end

function Base.show(io::IO,model::CompositeModel)
    print(io,string(typeof(model)),model.shape_model.components)
end

function volume_impl(model::CompositeModel,p,T,z,phase=:unknown,threaded=false,vol = vol0)
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

function saturation_pressure(model::CompositeModel,T::Real)
    if model.models.saturation isa SaturationModel
        method = SaturationCorrelation()
    else
        method = ChemPotVSaturation()
    end
    return saturation_pressure(model,T,method)
end

function saturation_pressure(cmodel::CompositeModel,T,method::SaturationMethod)
    model = cmodel.models
    nan = zero(T)/zero(T)
    psat,_,_ = saturation_pressure_impl(model.saturation,T,method)
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

export CompositeModel
