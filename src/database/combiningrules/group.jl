#implace base

function _group_sum!(out,groups::GroupParam,param)
    v = groups.n_groups_cache
    for (i,vi) in pairs(v)
        out[i] = dot(vi,param)
    end
    return out
end

function _group_sum!(out,groups::GroupParam,param::Number)
    v = groups.n_groups_cache
    for (i,vi) in pairs(v)
        out[i] = sum(vi)*param
    end
    return out
end

function group_sum!(out::SingleParameter,groups::GroupParam,param::SingleParameter)
    _group_sum!(out.values,groups,param)
    out.ismissingvalues .= param.ismissingvalues
    return out
end

function group_sum!(out,groups::GroupParam,param::Nothing)
    return _group_sum!(out,groups,true)
end

function group_sum!(out,groups::GroupParam,param)
    return _group_sum!(out,groups,param)
end

"""
    group_sum(groups::GroupParam,P::SingleParameter)

Given a `GroupParam` and a Single Parameter `P` for group data, it will return a single parameter `p` of component data, where:

pᵢ = ∑Pₖνᵢₖ

where `νᵢₖ` is the number of groups `k` at component `i`.

"""
function group_sum(groups::GroupParam,param::SingleParameter)
    gc = length(groups.components)
    out =  SingleParam(param.name,
                        groups.components,
                        zeros(float(eltype(param.values)),gc),
                        param.ismissingvalues,
                        param.sources,
                        param.sourcecsvs)
    return group_sum!(out,groups,param)
end

"""
    group_sum(groups::GroupParam,P::AbstractVector)

Given a `GroupParam` and a Vector `P` for group data, it will return a Vector `p` of component data, where:

pᵢ = ∑Pₖνᵢₖ

where `νᵢₖ` is the number of groups `k` at component `i`.

"""
function group_sum(groups::GroupParam,param::AbstractVector)
    out = similar(param,length(groups.components))
    return group_sum!(out,groups,param)
end

"""
    group_sum(groups::GroupParam,::Nothing)

Given a `GroupParam`, it will return a vector `p` of component data, where:

pᵢ = ∑νᵢₖ

where `νᵢₖ` is the number of groups `k` at component `i`.

"""
function group_sum(groups::GroupParam,::Nothing)
    out = zeros(Float64,length(groups.components))
    return group_sum!(out,groups,true)
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
    gc = length(groups.flattenedgroups)
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

function group_pairmean(f::T,groups::GroupParam,param::SingleOrPair) where {T}
    return SingleParam(param.name,groups.components,group_pairmean(f,groups,param.values))
end


function group_pairmean(f::T,groups::GroupParam,p::AbstractArray) where {T}
    _0 = zero(eltype(p))/one(eltype(p)) #we want a float type
    res = fill(_0,length(groups.components))
    return group_pairmean!(res,f,groups,p)
end

function group_pairmean!(res,f::T,groups::GroupParam,param::SingleOrPair) where {T}
    return group_pairmean!(res,f,groups,param.values)
end

function group_pairmean!(res,f,groups::GroupParam,p::AbstractMatrix)
    lgroups = 1:length(groups.flattenedgroups)
    lcomps = 1:length(res)
    _0 = zero(eltype(res))
    zz = groups.n_groups_cache
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

function group_pairmean!(res,f::T,groups::GroupParam,p::AbstractVector) where {T}
    lgroups = 1:length(groups.flattenedgroups)
    lcomps = 1:length(groups.components)
    zz = groups.n_groups_cache
    _0 = zero(eltype(res))
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
