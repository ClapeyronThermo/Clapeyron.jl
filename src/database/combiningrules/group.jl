"""
    group_sum(groups::GroupParam,P::SingleParameter)

Given a `GroupParam` and a Single Parameter `P` for group data, it will return a single parameter `p` of component data, where:

pᵢ = ∑Pₖνᵢₖ

where `νᵢₖ` is the number of groups `k` at component `i`.

"""
function group_sum(groups::GroupParam,param::SingleParameter)
    comp_vals = group_sum(groups,param.values)
    ismissingval = [ismissing(vi) for vi in comp_vals]
    return SingleParam(param.name,
                        groups.components,
                        comp_vals,
                        ismissingval,
                        param.sources,
                        param.sourcecsvs)
end

"""
    group_sum(groups::GroupParam,P::AbstractVector)

Given a `GroupParam` and a Vector `P` for group data, it will return a Vector `p` of component data, where:

pᵢ = ∑Pₖνᵢₖ

where `νᵢₖ` is the number of groups `k` at component `i`.

"""
function group_sum(groups::GroupParam,param::AbstractVector)
    v = groups.n_groups_cache
    return [dot(vi,param) for vi in v]
end

"""
    group_sum(groups::GroupParam,::Nothing)

Given a `GroupParam`, it will return a vector `p` of component data, where:

pᵢ = ∑νᵢₖ

where `νᵢₖ` is the number of groups `k` at component `i`.

"""
function group_sum(groups::GroupParam,::Nothing)
    v = groups.n_groups_cache
    return [sum(vi) for vi in v]
end

function group_sum(groups::GroupParam)
    return SingleParam("m",
                        groups.components,
                        group_sum(groups, nothing))
end
"""
    group_matrix(groups::GroupParam)

returns a matrix of size `(k,i)` with values νₖᵢ. when multiplied with a molar amount, it returns the amount of moles of each group.
"""
function group_matrix(groups::GroupParam)
    ng = groups.n_groups_cache
    comp = length(ng)
    gc = length(groups.i_flattenedgroups)
    return reshape(ng.v,(gc,comp))
end

"""
    group_pairmean(groups::GroupParam,param::PairParam)
    group_pairmean(f,groups::GroupParam,param::SingleParam)
Given a `GroupParam`and a parameter `P` it will return a single parameter `p` of component data, where:

pᵢ = ∑νᵢₖ(∑(νᵢₗ*P(i,j))) / ∑νᵢₖ(∑νᵢₗ)

where `νᵢₖ` is the number of groups `k` at component `i` and `P(i,j)` depends on the type of `P`:
- if `P` is a single paremeter, then `P(i,j) = f(P[i],P[j])`
- if `P` is a pair paremeter, then `P(i,j) = p[i,j]`

"""
function group_pairmean end

group_pairmean(groups::GroupParam,param) = group_pairmean(mix_mean,groups,param)

function group_pairmean(f::T,groups::GroupParam,param) where {T}
    return SingleParam(param.name,groups.components,group_pairmean(f,groups,param.values))
end

function group_pairmean(f,groups::GroupParam,p::AbstractMatrix)
    lgroups = 1:length(groups.i_flattenedgroups)
    lcomps = 1:length(groups.components)
    zz = groups.n_groups_cache
    _0 = zero(eltype(p))/one(eltype(p)) #we want a float type
    res = zeros(typeof(_0),length(lcomps))
    for i ∈ lcomps
        ẑ = zz[i]
        ∑ẑinv2 = 1/(sum(ẑ)^2)
        p_i = _0
        for k ∈ lgroups
            ẑk = ẑ[k]
            iszero(ẑk) && continue
            pk = p[k,k]
            p_i += ẑk*ẑk*pk
            for l ∈ 1:k - 1
                p_i += 2*ẑk*ẑ[l]*p[k,l]
            end
        end
        res[i] = p_i*∑ẑinv2
    end
    return res
end

function group_pairmean(f::T,groups::GroupParam,p::AbstractVector) where {T}
    lgroups = 1:length(groups.i_flattenedgroups)
    lcomps = 1:length(groups.components)
    zz = groups.n_groups_cache
    _0 = zero(eltype(p))
    res = zeros(eltype(p),length(lcomps))
    for i ∈ lcomps
        ẑ = zz[i]
        ∑ẑinv2 = 1/(sum(ẑ)^2)
        p_i = _0
        for k ∈ lgroups
            ẑk = ẑ[k]
            iszero(ẑk) && continue
            pk = p[k]
            p_i += ẑk*ẑk*pk
            for l ∈ 1:k - 1
                p_i += 2*ẑk*ẑ[l]*f(pk,p[l])
            end
        end
        res[i] = p_i*∑ẑinv2
    end
    return res
end


"""
    mix_segment!(groups::GroupParam,S = ones(length(@groups)),vst = ones(length(@groups)))

modifies implace the field `n_groups_cache` (`μᵢₖ`) in the `GroupParam`:
```
μᵢₖ =  νᵢₖ*Sₖ*vstₖ
```
Where `S` is a shape factor parameter for each group and `vst` is the segment size for each group.
used mainly for GC models (like `SAFTgammaMie`) in which the group fraction depends on segment size and shape factors.
"""
function mix_segment!(groups::GroupParam,s = ones(length(groups.flattenedgroups)),segment = ones(length(groups.flattenedgroups)))
    gc = 1:length(groups.flattenedgroups)
    comps = 1:length(groups.components)
    v = groups.n_groups_cache
    for k in 1:length(gc)
        for i in 1:length(comps)
            v[i][k] = v[i][k]*s[k]*segment[k]
        end
    end
    #SingleParam("mixed segment",groups.flattenedgroups,mixsegment,[false for i ∈ gc],String[],String[])
end
