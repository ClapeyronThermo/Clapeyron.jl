
include("GenericAncEvaluator.jl")
include("SaturationModel/SaturationModel.jl")
include("LiquidVolumeModel/LiquidVolumeModel.jl")
include("PolExpVapour.jl")
include("SolidModel/SolidHfus.jl")

Base.length(cmodel::CompositeModel) = length(cmodel.components)

function CompositeModel(components;
    liquid = RackettLiquid,
    gas = BasicIdeal,
    fluid=nothing,
    userlocations = String[],
    solid = nothing,
    saturation = LeeKeslerSat,
    melting = nothing,
    gas_userlocations = String[],
    liquid_userlocations = String[],
    fluid_userlocations = String[],
    solid_userlocations = String[],
    saturation_userlocations = String[],
    melting_userlocations = String[],
    verbose = false)

    if fluid !== nothing
        if fluid <: ActivityModel
            error("Activity models only represent the liquid phase. Please specify a gas phase model.")
        end
        gas = fluid
        liquid = fluid
        saturation = fluid
        gas_userlocations = fluid_userlocations
        liquid_userlocations = fluid_userlocations
        saturation_userlocations = fluid_userlocations
    end

    init_gas = init_model(gas,components,gas_userlocations,verbose)
    if typeof(liquid) <: EoSModel
        init_liquid = init_model(liquid,components,liquid_userlocations,verbose)
    else
        if liquid <: ActivityModel
            init_liquid = liquid(components;userlocations=liquid_userlocations,puremodel=gas,verbose)
        else
            init_liquid = init_model(liquid,components,liquid_userlocations,verbose)
        end
    end

    init_solid = init_model(solid,components,solid_userlocations,verbose)
    init_sat = init_model(saturation,components,saturation_userlocations,verbose)
    init_melt = init_model(melting,components,melting_userlocations,verbose)

    components = format_components(components)
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

function init_preferred_method(method::typeof(saturation_pressure),model::CompositeModel,kwargs)
    return init_preferred_method(saturation_pressure,model.saturation,kwargs)
end

function init_preferred_method(method::typeof(saturation_temperature),model::CompositeModel,kwargs)
    return init_preferred_method(saturation_temperature,model.saturation,kwargs)
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
    single_component_check(crit_pure,model)
    return crit_pure(model.saturation)
end

function x0_sat_pure(model::CompositeModel,T)
    p = x0_psat(model,T)
    vl = volume(model.liquid,p,T,phase=:l)
    vv = volume(model.gas,p,T,phase=:v)
    return vl,vv
end

function x0_psat(model::CompositeModel,T)
    ps,_,_ = saturation_pressure(model.saturation,T)
    return ps
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

#Michelsen TPFlash and rachford rice tpflash support
function init_preferred_method(method::typeof(tp_flash),model::CompositeModel,kwargs)
    return RRTPFlash(;kwargs...)
end

__tpflash_cache_model(model::CompositeModel,p,T,z) = PTFlashWrapper(model,T)

function PTFlashWrapper(model::CompositeModel,T::Number) 
    satmodels = split_model(model.saturation)
    gases = split_model(model.gas,1:length(model))
    sats = saturation_pressure.(satmodels,T)
    vv_pure = last.(sats)
    RT = R̄*T
    p_pure = first.(sats)
    μpure = only.(VT_chemical_potential_res.(gases,vv_pure,T))
    ϕpure = exp.(μpure ./ RT .- log.(p_pure .* vv_pure ./ RT))
    g_pure = [VT_gibbs_free_energy(gases[i],sats[i][2],T) for i in 1:length(model)]

    return PTFlashWrapper(model.components,model,sats,ϕpure,μpure)
end

function update_K!(lnK,wrapper::PTFlashWrapper{<:CompositeModel},p,T,x,y,volx,voly,phasex,phasey,β = nothing,inx = FillArrays.Fill(true,length(x)),iny = inx)
    model = wrapper.model
    sats = wrapper.sat
    #crits = wrapper.crit
    fug = wrapper.fug
    RT = R̄*T
    volx = volume(model.liquid, p, T, x, phase = phasex, vol0 = volx)
    lnϕy, voly = lnϕ(__gas_model(model), p, T, y; phase=phasey, vol0=voly)
    if is_vapour(phasey)
        for i in eachindex(lnK)
            if iny[i]
                ϕli = fug[i]
                p_i = sats[i][1]
                lnK[i] = log(p_i*ϕli/p) - lnϕy[i] + volx*(p - p_i)/RT
            end
        end
    else
        throw(error("CompositeModel does not support LLE equilibria."))
    end
    return lnK,volx,voly,NaN*one(T+p+first(x))
end

function __tpflash_gibbs_reduced(wrapper::PTFlashWrapper{<:CompositeModel},p,T,x,y,β,eq)
    return NaN*one(T+p+first(x))
end

function dgibbs_obj!(model::PTFlashWrapper{<:CompositeModel}, p, T, z, phasex, phasey,
    nx, ny, vcache, ny_var = nothing, in_equilibria = FillArrays.Fill(true,length(z)), non_inx = in_equilibria, non_iny = in_equilibria;
    F=nothing, G=nothing, H=nothing)
    throw(error("CompositeModel does not support gibbs energy optimization in MichelsenTPFlash."))
    #
end

function K0_lle_init(model::PTFlashWrapper{<:CompositeModel},p,T,z)
    throw(error("CompositeModel does not support LLE equilibria."))
end

export CompositeModel
