"""
    FluidCorrelation{V,L,Sat} <: RestrictedEquilibriaModel

Wrapper struct to signal that a `CompositeModel` uses correlations for calculation of saturation points, vapour and liquid phase volumes.
"""
struct FluidCorrelation{V,L,Sat} <: RestrictedEquilibriaModel
    components::Vector{String}
    gas::V
    liquid::L
    saturation::Sat
end

__gas_model(model::FluidCorrelation) = model.gas
activity_coefficient(model::FluidCorrelation, p, T,z=SA[1.]; phase = :unknown, threaded=true) = FillArrays.Ones(length(model)) 

function volume_impl(model::FluidCorrelation, p, T, z, phase, threaded, vol0)
    _0 = zero(p+T+first(z))
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

function crit_pure(model::FluidCorrelation)
    single_component_check(crit_pure,model)
    return crit_pure(model.saturation)
end

function x0_sat_pure(model::FluidCorrelation,T)
    p = x0_psat(model,T)
    vl = volume(model.liquid,p,T,phase=:l)
    vv = volume(model.gas,p,T,phase=:v)
    return vl,vv
end

function x0_psat(model::FluidCorrelation,T)
    ps,_,_ = saturation_pressure(model.saturation,T)
    return ps
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

__tpflash_cache_model(model::FluidCorrelation,p,T,z,equilibrium) = PTFlashWrapper(model,p,T,equilibrium)

function PTFlashWrapper(model::FluidCorrelation,p,T::Number,equilibrium::Symbol)
    satmodels = split_model(model.saturation)
    gases = split_model(model.gas,1:length(model))
    sats = saturation_pressure.(satmodels,T)
    vv_pure = last.(sats)
    RT = R̄*T
    p_pure = first.(sats)
    μpure = only.(VT_chemical_potential_res.(gases,vv_pure,T))
    ϕpure = exp.(μpure ./ RT .- log.(p_pure .* vv_pure ./ RT))
    g_pure = [VT_gibbs_free_energy(gases[i],sats[i][2],T) for i in 1:length(model)]
    return PTFlashWrapper(model.components,model,sats,ϕpure,μpure,equilibrium)
end

function update_K!(lnK,wrapper::PTFlashWrapper{<:FluidCorrelation},p,T,x,y,volx,voly,phasex,phasey,β = nothing,inx = FillArrays.Fill(true,length(x)),iny = inx)
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
        throw(error("Correlation-Based Composite Model does not support LLE equilibria."))
    end
    return lnK,volx,voly,NaN*one(T+p+first(x))
end

function __tpflash_gibbs_reduced(wrapper::PTFlashWrapper{<:FluidCorrelation},p,T,x,y,β,eq)
    return NaN*one(T+p+first(x))
end

function dgibbs_obj!(model::PTFlashWrapper{<:FluidCorrelation}, p, T, z, phasex, phasey,
    nx, ny, vcache, ny_var = nothing, in_equilibria = FillArrays.Fill(true,length(z)), non_inx = in_equilibria, non_iny = in_equilibria;
    F=nothing, G=nothing, H=nothing)
    throw(error("Correlation-Based Composite Model does not support gibbs energy optimization in MichelsenTPFlash."))
    #
end

function K0_lle_init(model::PTFlashWrapper{<:FluidCorrelation},p,T,z)
    throw(error("Correlation-Based Composite Model does not support LLE equilibria."))
end

function __eval_G_DETPFlash(model::PTFlashWrapper{<:FluidCorrelation},p,T,xi,equilibrium)
    throw(error("Correlation-Based Composite Model does not support DETPFlash."))
end

export FluidCorrelation
