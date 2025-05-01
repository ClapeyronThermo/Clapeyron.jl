#implace base
__get_group_sum_values(group::GroupParam) = group.n_flattenedgroups
__get_group_sum_values(group::MixedGCSegmentParam) = group.values

function _group_sum!(out,groups,param)
    v = __get_group_sum_values(groups)
    for (i,vi) in pairs(v)
        out[i] = dot(vi,param)
    end
    return out
end

function _group_sum!(out,groups,param::Number)
    v = __get_group_sum_values(groups)
    for (i,vi) in pairs(v)
        out[i] = sum(vi)*param
    end
    return out
end

function group_sum!(out::SingleParameter,groups,param::SingleParameter)
    _group_sum!(out.values,groups,param)
    v = __get_group_sum_values(groups)
    missingvals_comps = out.ismissingvalues
    missingvals_gc = param.ismissingvalues
    #173
    gc = length(v[1])
    comps = length(out.values)
    for i in 1:comps
        is_missing_i = false
        vi = v[i]
        for j in 1:gc
            if v[i] != 0 #if the group count is nonzero, then reduce the ismissing values
                is_missing_i = is_missing_i | missingvals_gc[j]
            end
        end
        missingvals_comps[i] = is_missing_i
    end
    return out
end

function group_sum!(out,groups,param::Nothing)
    return _group_sum!(out,groups,true)
end

function group_sum!(out,groups,param)
    return _group_sum!(out,groups,param)
end

"""
    group_sum(groups::GroupParam,P::SingleParameter)

Given a `GroupParam` and a Single Parameter `P` for group data, it will return a single parameter `p` of component data, where:

pᵢ = ∑Pₖνᵢₖ

where `νᵢₖ` is the number of groups `k` at component `i`.

"""
function group_sum(groups,param::SingleParameter{T}) where T
    gc = length(groups.components)
    out = SingleParam{T}(param.name,
                        groups.components,
                        zeros(T,gc),
                        fill(false,gc),
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
function group_sum(groups,param::AbstractVector)
    out = similar(param,length(groups.components))
    return group_sum!(out,groups,param)
end

"""
    group_sum(groups::GroupParam,::Nothing)

Given a `GroupParam`, it will return a vector `p` of component data, where:

pᵢ = ∑νᵢₖ

where `νᵢₖ` is the number of groups `k` at component `i`.

"""
function group_sum(groups,::Nothing)
    v = __get_group_sum_values(groups)
    out = zeros(eltype(v[1]),length(groups.components))
    return group_sum!(out,groups,true)
end

function group_sum(groups)
    return SingleParam("m",
                        groups.components,
                        group_sum(groups, nothing))
end
"""
    group_matrix(groups::MixedGCSegmentParam)

returns a matrix of size `(k,i)` with values νₖᵢ. when multiplied with a molar amount, it returns the amount of moles of each group.
"""
function group_matrix(groups::MixedGCSegmentParam)
    vals = groups.values
    comp = length(vals)
    gc = length(vals[1])
    return reshape(vals.v,(gc,comp))
end

function group_fractions(groups::MixedGCSegmentParam{T1},z::AbstractVector{T2}) where {T1,T2}
    v = groups.values
    ng = length(v[1])
    x = similar(z, Base.promote_type(T1,T2), ng)
    M = group_matrix(groups)
    mul!(x,M,z)
    return x
end

function group_fractions(groups::GroupParam,z)
    ng = length(groups.flattenedgroups)
    n_flattenedgroups = groups.n_flattenedgroups
    x = similar(z,ng)
    nc = length(z)
    fill!(x,zero(eltype(x)))
    @inbounds for i in 1:nc
        n_i = n_flattenedgroups[i]
        z_i = z[i]
        for k in 1:ng
            x[k] += z_i*n_i[k]
        end
    end
    return x
end

"""
    group_pairmean(groups::GroupParam,param::PairParam)
    group_pairmean(f,groups::GroupParam,param::SingleParam)
Given a `GroupParam` and a parameter `P` it will return a single parameter `p` of component data, where:

pᵢ = ∑νᵢₖ(∑(νᵢₗ*P(i,j))) / ∑νᵢₖ(∑νᵢₗ)

where `νᵢₖ` is the number of groups `k` at component `i` and `P(i,j)` depends on the type of `P`:
- if `P` is a single paremeter, then `P(i,j) = f(P[i],P[j])`
- if `P` is a pair paremeter, then `P(i,j) = p[i,j]`

"""
function group_pairmean end

group_pairmean(groups,param) = group_pairmean(mix_mean,groups,param)

function group_pairmean(f::T,groups,param::SingleOrPair) where {T}
    return SingleParam(param.name,groups.components,group_pairmean(f,groups,param.values))
end

function group_pairmean(f::F,groups,p::AbstractArray) where {F}
    v = __get_group_sum_values(groups)
    T = Base.promote_eltype(1.0,p,v[1])
    res = zeros(T, length(groups.components))
    return group_pairmean!(res,f,groups,p)
end

function group_pairmean!(res,f::F,groups,param::SingleOrPair) where {F}
    return group_pairmean!(res,f,groups,param.values)
end

function group_pairmean!(res,f,groups,p::AbstractMatrix)
    zz = __get_group_sum_values(groups)
    lgroups = 1:length(zz[1])
    lcomps = 1:length(res)
    _0 = zero(eltype(res))
    
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

function group_pairmean!(res,f::T,groups,p::AbstractVector) where {T}
    zz = __get_group_sum_values(group)
    lgroups = 1:length(zz[1])
    lcomps = 1:length(groups.components)
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
    mix_segment!(groups::MixedGCSegmentParam,S = ones(length(@groups)),vst = ones(length(@groups)))

modifies implace the field `n_groups_cache` (`μᵢₖ`) in the `GroupParam`:
```
μᵢₖ = νᵢₖ*Sₖ*vstₖ
```
Where `S` is a shape factor parameter for each group and `vst` is the segment size for each group.
used mainly for GC models (like `SAFTgammaMie`) in which the group fraction depends on segment size and shape factors.
"""
function mix_segment!(groups,s = ones(length(groups.flattenedgroups)),segment = ones(length(groups.flattenedgroups)))
    v = __get_group_sum_values(groups)
    ng = length(v[1])
    nc = length(v)
    for i in 1:nc
        vi = v[i]
        for k in 1:ng
            vi[k] = vi[k]*s[k]*segment[k]
        end
    end
    return groups
    #SingleParam("mixed segment",groups.flattenedgroups,mixsegment,[false for i ∈ gc],String[],String[])
end


function group_pairmean2(groups::GroupParameter,param::PairParam)
    newvals = group_pairmean2!(groups,copy(param.values))
    return PairParam(param.name,groups.components,newvals,fill(false,size(newvals)),param.sources,param.sourcecsvs)
end

function group_pairmean2!(groups,mat)
    l_gc = length(groups.flattenedgroups)
    l_c = length(groups.components)
    _0 = zero(eltype(mat))
    newmat = fill(_0,(l_c,l_c))
    n = groups.n_flattenedgroups
    for i ∈ 1:l_c
        for j ∈ 1:l_c
            res = _0
            sumn = _0
            for k in 1:l_gc
                for l in 1:l_gc
                    res += n[i][k]*n[j][l]*mat[k,l]
                    sumn += n[i][k]*n[j][l]
                end
            end
            newmat[i,j] = res/sumn
        end
    end
    return newmat
end
