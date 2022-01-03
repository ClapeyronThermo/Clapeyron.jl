


struct CachedEoS{EoS<:EoSModel,S} <: EoSModel
    model::EoS
    cache::Dict{Symbol,Any}
    saturation_pressure_aprox::S
end

const CACHED_SAT_PURE_APROX{S} =  CachedEoS{<:Any,S}


#is 256 too low? to high?
sizehint(::CACHED_SAT_PURE_APROX{Nothing}) = 256

function saturation_pressure_aprox_strategy(model::CachedEoS)
    return model.saturation_pressure_aprox
end
Base.length(model::CachedEoS) = Base.length(model.model)

function CachedEoS(model::EoSModel;saturation_pressure_aprox = nothing)
    cache=Dict{Symbol,Any}()
    return CachedEoS(model,cache,saturation_pressure_aprox)
end

function eos(model::CachedEoS,V,T,z=SA[1.0],phase=:unknown)
    return eos(model.model,V0,T0)
end

function eos_res(model::CachedEoS,V,T,z=SA[1.0],phase=:unknown)
    return eos_res(model.model,V,T)
end
mw(model::CachedEoS) = mw(model.model)

function Base.show(io::IO,mime::MIME"text/plain",model::CachedEoS)
    println(io,"EoS Model with cached results")
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

function split_model(model::CachedEoS)
    if haskey(model.cache,:split_model)
        return model.cache[:split_model]
    else
        primal_splitted_model = split_model(model.model)
        res = [CachedEoS(modeli,saturation_pressure_aprox = model.saturation_pressure_aprox) for modeli in primal_splitted_model]
        model.cache[:split_model] = res
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
function saturation_pressure(model::CachedEoS{<:Any},T::Real)
    aprox_strategy = saturation_pressure_aprox_strategy(model)
    return saturation_pressure_cached(aprox_strategy,model,T)
end

#no caching
function saturation_pressure_cached(::Nothing,model::CachedEoS,T)
    return saturation_pressure(model.model,T)
end


"""
    saturation_pressure!(model::CachedEoS,T,x0=x0_sat_pure(model,T))
    performs a calculation of [`saturation_pressure`](@ref) and stores the results in a vector, sorted by temperature.

"""
function saturation_pressure!(model::CachedEoS,T::S,x0 = x0_sat_pure(model,T)) where S<:Real
    #empty cache, create one
    if !haskey(model.cache,:saturation_pressure)
        (T_c, p_c, V_c) = crit_pure(model)  
        T_TYPE = promote_type(S,typeof(T_c))
        T_vec,P_vec,Vl_vec,Vv_vec = T_TYPE[],T_TYPE[],T_TYPE[],T_TYPE[]
        hint = sizehint(model) #gives an expected length
        sizehint!(T_vec,hint)
        sizehint!(P_vec,hint)
        sizehint!(Vl_vec,hint)
        sizehint!(Vv_vec,hint)
        push!(T_vec,T_c)
        push!(P_vec,p_c)
        push!(Vl_vec,V_c)
        push!(Vv_vec,V_c)
        model.cache[:saturation_pressure] = (T_vec,P_vec,Vl_vec,Vv_vec)
    end
    T_vec,P_vec,Vl_vec,Vv_vec = model.cache[:saturation_pressure]
    position = searchsorted(T_vec,T)
    #one and only one match, that means that T was already evaluated
    if first(position) == last(position)
        idx = first(position)
        _P0 = convert(typeof(T),P_vec[idx])
        _Vl = convert(typeof(T),Vl_vec[idx])
        _Vv = convert(typeof(T),Vv_vec[idx])
        return _P0,_Vl,_Vv
    end
    #calculate sat pure
    T_c = last(T_vec)  
    if T_c < T
        insert_first = true
    end
    P0,Vl,Vv = saturation_pressure(model.model,T,x0_sat_pure(model,T))
    #we dont want to store points too close
    insert_values = true
    if length(T_vec) >= sizehint(model)
        insert_values = false
    end

    #avoid splitting if possible
    insert_first = false
    Tmin = first(T_vec)  
    if Tmin > T
        insert_first = true
    end
    if insert_values
        if insert_first
        pushfirst!(T_vec,T)
        pushfirst!(P_vec,P0)
        pushfirst!(Vl_vec,Vl)
        pushfirst!(Vv_vec,Vv)
        else
        splice!(T_vec, position, [T]) 
        splice!(P_vec, position, [P0])     
        splice!(Vl_vec, position, [Vl]) 
        splice!(Vv_vec, position, [Vv])  
        end   
    end
    return P0,Vl,Vv
end
include("sat_pure_aprox/TPInterpolation.jl")

export CachedEoS
