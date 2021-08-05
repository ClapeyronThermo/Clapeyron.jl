
abstract type AbstractTPInterpolation <: SatPureAproximation
struct invTlogPInterpolation 
    max_length::Int
end

function sizehint(cachedmodel::CACHED_SAT_PURE_APROX{<:AbstractTPInterpolation}) = 256
    return cachedmodel.sat_pure_aprox.max_length
end

const TPInterpolation = invTlogPInterpolation

function TPInterpolation() = TPInterpolation(256) #max length of 200 data points
end

function sat_pure_cached(::TPInterpolation,model::CachedEoS,T)
    res = sat_pure(model.model,T)
end


