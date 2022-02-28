struct UNIFACCache <: EoSModel
    components::Vector{String}
    r::Vector{Float64}
    q::Vector{Float64}
    q_p::Vector{Float64}
    m::Vector{Float64}
end

UNIFACCache(groups::GroupParam,params) = UNIFACCache(groups,params.Q.values,params.R.values)

function UNIFACCache(groups::GroupParam,Q,R)
    comps = 1:length(groups.components)
    v = groups.n_flattenedgroups
    r = [dot(vi,R) for vi in v]
    q = [dot(vi,Q) for vi in v]
    q_p = r.^(3/4)
    comp_segment = zeros(length(comps))
    for i ∈ comps
        res_i = 0.0
        vi = v[i]
        groups_i = groups.i_groups[i]
        for idx ∈ 1:length(groups_i)
            k = groups_i[idx]
            res_i += vi[k]
        end
        comp_segment[i] = res_i
    end
    m = comp_segment
    return UNIFACCache(groups.components,r,q,q_p,m)
end



