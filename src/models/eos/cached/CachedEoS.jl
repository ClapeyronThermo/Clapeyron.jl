struct CachedEoS{EoS<:EoSModel} <: EoSModel
    model::EoS
    cache::Dict{Symbol,Any}
end

Base.length(model::CachedEoS) = Base.length(model.model)

function CachedEoS(model::EoSModel)
    cache=Dict{Symbol,Any}()
    return CachedEoS(model,cache)
end

function eos(model::CachedEoS,V,T,z=SA[1.0],phase=:unknown)
    return eos(model.model,V0,T0)
end

function eos_res(model::CachedEoS,V,T,z=SA[1.0],phase=:unknown)
    return eos_res(model.model,V,T)
end
mw(model::CachedEoS) = mw(model.model)

function Base.show(io::IO,mime::MIME"text/plain",model::CachedEoS)
    println(io,"EoS Model with cached results:")
    show(io,mime,model.model)
end

function Base.show(io::IO,model::CachedEoS)
    print(io,"CachedEoS(",model.model,")")
end

function lb_volume(model::CachedEoS,z=SA[1.0])
    return lb_volume(model.model)
end

function T_scale(model::CachedEoS,z=SA[1.0])
    return T_scale(model.model)
end

function p_scale(model::CachedEoS,z=SA[1.0])
    return p_scale(model.mode)
end

function crit_pure(model::CachedEoS)
    if haskey(model.cache,:crit_pure)
        return model.cache[:crit_pure]
    else
        local res
        try
            res = crit_pure(model.model)
        catch
            res = (NaN,NaN,NaN)
        end
        model.cache[:crit_pure] = res
        return res
    end
end

function acentric_factor(model::CachedEoS)
    if haskey(model.cache,:acentric_factor)
        return model.cache[:acentric_factor]
    else
        local res
        try
            res = acentric_factor(model.model)
        catch
            res = NaN
        end
        model.cache[:acentric_factor] = res
        return res
    end
end

#for some reason, it does not work without an overloading here
function sat_pure(model::CachedEoS,T::Real)
   return sat_pure(model.model,T)
end

export CachedEoS
