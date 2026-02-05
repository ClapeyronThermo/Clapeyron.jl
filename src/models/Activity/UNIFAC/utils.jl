#a GC averaged UNIFAC.
struct UNIFACCache{T} <: EoSModel
    components::Vector{String}
    r::Vector{T}
    q::Vector{T}
    q_p::Vector{T}
end

UNIFACCache(components,r,q,q_p) = UNIFACCache{eltype(q_p)}(components,r,q,q_p)

UNIFACCache(groups::GroupParam,params) = UNIFACCache(groups,params.Q,params.R)

function UNIFACCache(groups::GroupParam,Q,R)
    r = group_sum(groups,R.values)
    q = group_sum(groups,Q.values)
    q_p = r.^(3/4)
    return UNIFACCache(groups.components,r,q,q_p)
end

function recombine_unifac_cache!(cache::UNIFACCache,groups,params)
    Q = params.Q
    R = params.R
    group_sum!(cache.r,groups,R.values)
    group_sum!(cache.q,groups,Q.values)
    cache.q_p .= cache.r.^(3/4)
    return cache
end
