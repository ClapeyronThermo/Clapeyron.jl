"""
    FluidCorrelation{V,L,Sat,Cp} <: RestrictedEquilibriaModel

Wrapper struct to signal that a `CompositeModel` uses correlations for calculation of saturation points, vapour and liquid phase volumes.
"""
struct FluidCorrelation{V,L,Sat,Cp} <: RestrictedEquilibriaModel
    components::Vector{String}
    gas::V
    liquid::L
    saturation::Sat
    liquid_cp::Cp
end

function FluidCorrelation(_components;
                            gas_volume = nothing,
                            liquid_volume = nothing,
                            saturation = nothing,
                            liquid_cp = nothing
                            ,gas_volume_userlocations = String[],
                            liquid_volume_userlocations = String[],
                            saturation_userlocations = String[],
                            liquid_cp_userlocations = String[],
                            verbose = false,
                            liquid_reference_state = :ntp) #=Coolprop uses this reference=#

    components = format_components(_components)
    if gas_volume !== nothing
        init_gas = init_model(gas_volume,components,gas_volume_userlocations,verbose)
    else
        init_gas = nothing
    end

    if liquid_volume !== nothing
        init_liquid = init_model(liquid_volume,components,liquid_volume_userlocations,verbose)
    else
        init_liquid = nothing
    end

    if saturation !== nothing
        init_sat = init_model(saturation,components,saturation_userlocations,verbose)
    else
        init_sat = nothing
    end

    if liquid_cp !== nothing
        init_cp = init_model(liquid_cp,components,liquid_cp_userlocations,verbose)
    else
        init_cp = nothing
    end


    model = FluidCorrelation(components,init_gas,init_liquid,init_sat,init_cp)
    return model
end

function Base.show(io::IO,mime::MIME"text/plain",model::FluidCorrelation)
    print(io,"Fluid Correlation Model")
    length(model) == 1 && print(io, " with 1 component:")
    length(model) > 1 && print(io, " with ", length(model), " components:")
    model.gas !== nothing && print(io,'\n',"Gas Model: ",model.gas)
    model.liquid !== nothing && print(io,'\n',"Liquid Model: ",model.liquid)
    model.saturation !== nothing && print(io,'\n',"Saturation Model: ",model.saturation)
    model.liquid_cp !== nothing && print(io,'\n',"Liquid Caloric Model: ",model.liquid_cp)
end

reference_state(model::FluidCorrelation) = reference_state(model.gas)
Base.eltype(model::FluidCorrelation) = Base.promote_eltype(__γ_unwrap(model),gas_model(model))

function idealmodel(model::FluidCorrelation{V}) where V
    idealmodel(model.gas)
end

@inline gas_model(Base.@specialize(model::FluidCorrelation)) = model.gas
liquid_model(model::FluidCorrelation) = model.liquid

function PT_property(model::FluidCorrelation,p,T,z,phase,threaded,vol0,f::F,vol::V) where {F,V}
    if is_vapour(phase)
        return PT_property(model.gas,p,T,z,phase,threaded,vol0,f,vol)
    else #liquid or unknown
        wrapper = PTFlashWrapper(model,p,T,z,:vle)
        return PT_property(wrapper,p,T,z,phase,threaded,vol0,vol)
    end
end

struct FluidCorrelationSaturation{M} <: SaturationMethod
    method::M
end

function init_preferred_method(method::typeof(saturation_pressure),model::FluidCorrelation,kwargs)
    inner = init_preferred_method(method,model.saturation,kwargs)
    return FluidCorrelationSaturation(inner)
end

function init_preferred_method(method::typeof(saturation_temperature),model::FluidCorrelation,kwargs)
    inner = init_preferred_method(method,model.saturation,kwargs)
    return FluidCorrelationSaturation(inner)
end

function saturation_pressure_ad2(result,model::M,T::ForwardDiff.Dual) where M <: FluidCorrelation
    p = if has_a_res(model.saturation) #using an EoSModel as saturation provider
        #NOTE: we are using the volumes of the correlations instead of the volumes calculated via saturation_pressure of the helmholtz saturation model.
        #i don't know how correct is this, if you find an error on derivatives while calculating bulk properties of composite models and read this,
        #open a issue.
        first(saturation_pressure_ad2(result,model.saturation,T))
    else
        first(saturation_pressure(model.saturation,T)) #AD though method, directly
    end

    liq = model.liquid
    z = SA[1.0]
    vl = if has_a_res(liq)
        tup = (liq,p,T,z)
        λtup = (liq,primalval(p),primalval(T),z)
        volume_ad(result[2],tup,λtup)
    else
        volume(liq,p,T,z,phase = :l)
    end

    gas = model.gas
    vv = if is_idealmodel(gas)
        Rgas(gas)*T/p
    else
        tup = (gas,p,T,z)
        λtup = (gas,primalval(p),primalval(T),z)
        volume_ad(result[3],tup,λtup)
    end
    return p,vl,vv
end

function saturation_pressure_ad(__model::M,result::RES,tup::TUP1,tup_primal::TUP2) where {M <: FluidCorrelation,RES,TUP1,TUP2}
    dm,∂T = tup
    m0,λT = tup_primal
    λp,λvl,λvv = result
    if any(has_dual,tup) # do check here to avoid recomputation of pressure if no AD
        if has_a_res(__model.saturation) #using an EoSModel as saturation provider
            vl0 = volume(m0.saturation,λp,λT,phase = :l)
            vv0 = volume(m0.saturation,λp,λT,phase = :v)
            ∂p = saturation_pressure_ad(__model.saturation,(λp,vl0,vv0),(dm.saturation,∂T),(m0.saturation,λT))
        else
            ∂p = first(saturation_pressure(dm.saturation,λT)) #AD though method, directly
        end

        ∂liq = dm.liquid
        λliq = m0.liquid
        z = SA[1.0]
        ∂vl = if has_a_res(λliq)
            tup_l = (λliq,∂p,∂T,z)
            λtup_l = (∂liq,λp,λT,z)
            volume_ad(λvl,tup_l,λtup_l)
        else
            volume(∂liq,∂p,∂T,z,phase = :l)
        end

        ∂gas = dm.gas
        λgas = m0.gas
        ∂vv = if is_idealmodel(λgas)
            Rgas(λgas)*∂T/∂p
        else
            tup_v = (∂gas,∂p,∂T,z)
            λtup_v = (gas,λp,λT,z)
            volume_ad(λvv,tup_v,λtup_v)
        end
    
        return ∂p,∂vl,∂vv
    end
    return result
end

function saturation_temperature_ad(__model::M,result,tup,tup_primal) where M <: FluidCorrelation
    dm,∂p = tup
    m0,λp = tup_primal
    λT,λvl,λvv = result
    if any(has_dual,tup) # do check here to avoid recomputation of pressure if no AD
        if has_a_res(__model.saturation) #using an EoSModel as saturation provider
            vl0 = volume(m0.saturation,λp,λT,phase = :l)
            vv0 = volume(m0.saturation,λp,λT,phase = :v)
            ∂p = saturation_temperature_ad(__model.saturation,(λT,vl0,vv0),(dm.saturation,∂p),(m0.saturation,λp))
        else
            ∂p = saturation_temperature_corr_ad(λp,(dm.saturation,∂p),(m0.saturation,λp))
        end

        ∂liq = dm.liquid
        λliq = m0.liquid
        z = SA[1.0]
        ∂vl = if has_a_res(λliq)
            tup_l = (λliq,∂p,∂T,z)
            λtup_l = (∂liq,λp,λT,z)
            volume_ad(λvl,tup_l,λtup_l)
        else
            volume(∂liq,∂p,∂T,z,phase = :l)
        end

        ∂gas = dm.gas
        λgas = m0.gas
        ∂vv = if is_idealmodel(λgas)
            Rgas(λgas)*∂T/∂p
        else
            tup_v = (∂gas,∂p,∂T,z)
            λtup_v = (gas,λp,λT,z)
            volume_ad(λvv,tup_v,λtup_v)
        end
    
        return ∂T,∂vl,∂vv
    end
    return result
end

__γ_unwrap(model::FluidCorrelation) = IdealLiquidSolution()

reference_chemical_potential_type(model::FluidCorrelation) = :zero

function volume_impl(model::FluidCorrelation, p, T, z, phase, threaded, vol0)
    _0 = zero(p+T+first(z))
    _1 = one(_0)
    if model.gas === nothing && model.liquid !== nothing
        return _1*volume_impl(model.liquid,p,T,z,phase,threaded,vol0)
    elseif model.liquid === nothing && model.gas !== nothing
        return _1*volume_impl(model.gas,p,T,z,phase,threaded,vol0)
    end

    #nan = _0/_0
    if is_liquid(phase)
        return volume(model.liquid, p, T, z; phase, threaded, vol0)
    elseif is_vapour(phase)
        return volume(model.gas, p, T, z; phase, threaded, vol0)
    else
        wrapper = PTFlashWrapper(model,p,T,z,:vle)
        return volume(wrapper,p,T,z;phase,threaded,vol0)
    end
end

function saturation_pressure_impl(model::FluidCorrelation,T,methodwrapper::FluidCorrelationSaturation{M}) where M
    method = methodwrapper.method
    nan = zero(T)/zero(T)
    psat,_,_ = saturation_pressure(model.saturation,T,method)
    if !isnan(psat)
        vl = volume(model.liquid,psat,T,phase=:l)
        vv = volume(model.gas,psat,T,phase=:v)
        sat = psat,vl,vv
    #if psat fails, there are two options:
    #1- over critical point -> nan nan nan
    #2- saturation failed -> nan nan nan
    else
        sat = nan,nan,nan
    end
    return sat
end

function crit_pure(model::FluidCorrelation)
    single_component_check(crit_pure,model)
    return crit_pure(model.saturation)
end

function saturation_temperature_impl(model::FluidCorrelation,p,methodwrapper::FluidCorrelationSaturation{M}) where M
    method = methodwrapper.method
    nan = zero(p)/zero(p)
    Tsat,_,_ = saturation_temperature(model.saturation,p,method)
    if !isnan(Tsat)
        vl = volume(model.liquid,p,Tsat,phase=:l)
        vv = volume(model.gas,p,Tsat,phase=:v)
        sat = Tsat,vl,vv
    #if psat fails, there are two options:
    #1- over critical point -> nan nan nan
    #2- saturation failed -> nan nan nan
    else
        sat = nan,nan,nan
    end
    return sat
end

function dpdT_saturation(model::FluidCorrelation,v1::Number,v2,T)
    return dpdT_saturation(model.saturation,v1,v2,T)
end

function init_preferred_method(method::typeof(tp_flash),model::FluidCorrelation,kwargs)
    RRTPFlash(;kwargs...)
end

function init_preferred_method(method::typeof(tp_flash),model::FluidCorrelation{<:IdealModel},kwargs)
    RRTPFlash(;nacc = 0,kwargs...)
end

function __tpflash_cache_model(model::FluidCorrelation,p,T,z,equilibrium)
    PTFlashWrapper(model,p,T,z,equilibrium)
end

function ∂lnϕ_cache(model::PTFlashWrapper{FluidCorrelation{<:IdealModel}}, p, T, z, dt::Val{B}) where B
    return nothing
end

export FluidCorrelation
