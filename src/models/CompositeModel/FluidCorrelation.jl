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
    if gas_volume != nothing
        init_gas = init_model(gas_volume,components,gas_volume_userlocations,verbose)
    else
        init_gas = nothing
    end

    if liquid_volume != nothing
        init_liquid = init_model(liquid_volume,components,liquid_volume_userlocations,verbose)
    else
        init_liquid = nothing
    end

    if saturation != nothing
        init_sat = init_model(saturation,components,saturation_userlocations,verbose)
    else
        init_sat = nothing
    end

    if liquid_cp != nothing
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

function idealmodel(model::FluidCorrelation{V}) where V
    idealmodel(model.gas)
end

gas_model(model::FluidCorrelation) = model.gas

function PT_property(model::FluidCorrelation,p,T,z,phase,threaded,vol0,f::F,USEP::Val{UseP}) where {F,UseP}
    if is_unknown(phase) || phase == :stable
        #only liquid phase available
        if model.liquid_cp === nothing && model.gas !== nothing
            return PT_property(gas_model(model),p,T,z,phase,threaded,vol0,f,USEP)
        end

        #only gas phase available
        if model.gas === nothing && model.liquid_cp !== nothing
            return PT_property(model.liquid_cp,p,T,z,phase,threaded,vol0,f,USEP)
        end

        #single component, both phase models available, use saturation pressure to select
        if model.gas !== nothing && model.liquid_cp !== nothing && length(model) == 1
            psat,_,_ = saturation_pressure(model.saturation,T)
            if isnan(psat)
                #fluid correlations are only available on subcritical regime, where
                #we can distinguish a liquid from a vapour.
                nan = zero(Base.promote_eltype(model,p,T,z))
                return PT_property(model.gas,nan,nan,z,phase,threaded,vol0,f,USEP)
            end

            #phase identified, call again with correct phase
            new_phase = p > psat ? :l : :v
            return PT_property(model,p,T,z,new_phase,threaded,vol0,f,USEP)
        end
        throw(error("multicomponent automatic phase detection not implemented for $(typeof(model))"))
    end

    if is_liquid(phase)
        #for bulk properties that arent volume, the liquid_cp model contains a valid helmholtz model
        return PT_property(model.liquid_cp,p,T,z,phase,threaded,vol0,f,USEP)
    elseif is_vapour(phase)
        return PT_property(model.gas,p,T,z,phase,threaded,vol0,f,USEP)
    else
        throw(error("invalid phase specifier for FluidCorrelation: $phase"))
    end
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

    nan = _0/_0
    if is_liquid(phase)
        return volume(model.liquid, p, T, z; phase, threaded, vol0)
    elseif is_vapour(phase)
        return volume(model.gas, p, T, z; phase, threaded, vol0)
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
                        return volume(model.liquid, p, T, z; phase, threaded, vol0)
                    else #gas phase
                        return volume(model.gas, p, T, z; phase, threaded, vol0)
                    end
                else #something failed on saturation_pressure, not related to passing the critical point
                    return nan
                end

            end

            return nan
        else
            throw(error("multicomponent automatic phase detection not implemented for $(typeof(model))"))
            return nan
        end
    end
end

function init_preferred_method(method::typeof(saturation_pressure),model::FluidCorrelation,kwargs)
    return init_preferred_method(method,model.saturation,kwargs)
end

function init_preferred_method(method::typeof(saturation_temperature),model::FluidCorrelation,kwargs)
    return init_preferred_method(method,model.saturation,kwargs)
end

function saturation_pressure(model::FluidCorrelation,T,method::SaturationMethod)
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

function x0_sat_pure(model::FluidCorrelation,T,crit = nothing)
    p = x0_psat(model,T,crit)
    vl = volume(model.liquid,p,T,phase=:l)
    vv = volume(model.gas,p,T,phase=:v)
    return vl,vv
end

function x0_psat(model::FluidCorrelation,T,crit = nothing)
    ps,_,_ = saturation_pressure(model.saturation,T)
    return ps
end

function crit_pure(model::FluidCorrelation)
    single_component_check(crit_pure,model)
    return crit_pure(model.saturation)
end

function saturation_temperature(model::FluidCorrelation,p,method::SaturationMethod)
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

function init_preferred_method(method::typeof(tp_flash),model::FluidCorrelation,kwargs)
    RRTPFlash(;kwargs...)
end

function init_preferred_method(method::typeof(tp_flash),model::FluidCorrelation{<:IdealModel},kwargs)
    RRTPFlash(;nacc = 0,kwargs...)
end

function __tpflash_cache_model(model::FluidCorrelation,p,T,z,equilibrium) 
    PTFlashWrapper(model,p,T,equilibrium)
end

function ∂lnϕ_cache(model::PTFlashWrapper{FluidCorrelation{<:IdealModel}}, p, T, z, dt::Val{B}) where B
    return nothing
end

function __tpflash_gibbs_reduced(wrapper::PTFlashWrapper{<:FluidCorrelation},p,T,x,y,β,eq,vols)
    return NaN*one(T+p+first(x))
end

function K0_lle_init(model::PTFlashWrapper{<:FluidCorrelation},p,T,z)
    throw(error("Correlation-Based Composite Model does not support LLE equilibria."))
end

function __eval_G_DETPFlash(model::PTFlashWrapper{<:FluidCorrelation},p,T,xi,equilibrium)
    throw(error("Correlation-Based Composite Model does not support DETPFlash."))
end

export FluidCorrelation
