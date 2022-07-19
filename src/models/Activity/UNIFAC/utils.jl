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



