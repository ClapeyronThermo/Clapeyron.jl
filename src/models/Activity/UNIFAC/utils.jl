#a GC averaged UNIFAC.
struct UNIFACCache <: EoSModel
    components::Vector{String}
    r::Vector{Float64}
    q::Vector{Float64}
    q_p::Vector{Float64}
    m::Vector{Float64}
end

UNIFACCache(groups::GroupParam,params) = UNIFACCache(groups,params.Q,params.R)

function UNIFACCache(groups::GroupParam,Q,R)
    r = group_sum(groups,R.values)
    q = group_sum(groups,Q.values)
    q_p = r.^(3/4)
    m = group_sum(groups,nothing)
    return UNIFACCache(groups.components,r,q,q_p,m)
end

function recombine_unifac_cache!(cache::UNIFACCache,groups,params)
    Q = params.Q
    R = params.R
    group_sum!(cache.r,groups,R.values)
    group_sum!(cache.q,groups,Q.values)
    cache.q_p .= cache.r.^(3/4)
    group_sum!(cache.m,groups,nothing)
    return cache
end
