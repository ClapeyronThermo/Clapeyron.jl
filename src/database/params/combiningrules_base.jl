"""
    kij_mix(f,p::ClapeyronParam,k::PairParam)::PairParam
    kij_mix(f,p::ClapeyronParam)::PairParam

General combining rule for pair parameter with a `kᵢⱼ` interaction parameter. returns a pair parameter with non diagonal entries equal to:
```
pᵢⱼ = f(pᵢ,pⱼ,kᵢⱼ)
```
Where `f` is a 'combining' function that follows the rules:
```
f(pᵢ,pⱼ,0) = f(pⱼ,pᵢ,0)
f(pᵢ,pᵢ,0) = pᵢ
```
and `k` must follow:
```
kᵢᵢ = 0 
```
Ignores non-diagonal entries already set.

If a Single Parameter is passed as input, it will be converted to a Pair Parameter with `pᵢᵢ = pᵢ`.
"""
function kij_mix(f::F,param::SingleOrPair) where F
    out = PairParam(param) #copy the input
    return kij_mix!(f,out,nothing)
end

function kij_mix!(f::F,out::PairParameter) where F
    return kij_mix!(f,out,nothing)
end

function kij_mix!(f::F,out::PairParameter,::Nothing) where F
    N = length(out.components)
    k = FillArrays.Zeros(N,N)
    out_missing = out.ismissingvalues
    kij_mix!(f,out.values,k,out_missing)
    #if kij is missing, then the output values should be the same as the input values.
    #no missing prop has to be done
    return out
end

function kij_mix!(f::F,out::PairParameter,K::PairParameter) where F
    out_missing = out.ismissingvalues
    kij_mix!(f,out.values,K.values,out_missing)
    #should consider the two.
    out_missing .= out_missing .& K.ismissingvalues
    #but diagonals are all non-missing, by default:
    for i in diagind(out_missing)
        out_missing[i] = false
    end
    return out
end

#p,K: matrices
#B: bool matrix
function kij_mix!(f::F,p,K,B) where F
    N = LinearAlgebra.checksquare(p)
    for j ∈ 1:N
        p_j = p[j,j]
        for i ∈ 1:N
            if B[j,i]
                p_i = p[i,i]
                p_ji = f(p_i,p_j,K[j,i])
                p[j,i] = p_ji
            end
        end
    end
end



"""
    pair_mix(g,P::ClapeyronParam,Q::ClapeyronParam)::PairParam
    pair_mix(g,P::ClapeyronParam,Q::ClapeyronParam)::PairParam

General combining rule for a pair and a single parameter. returns a pair parameter `P` with non diagonal entries equal to:
```
Pᵢⱼ = g(Pᵢ,Pⱼ,Qᵢ,Qⱼ,Qᵢⱼ)
```
Where `f` is a 'combining' function that follows the rules:
```
Pᵢⱼ = Pⱼᵢ = g(Pᵢ,Pⱼ,Qᵢ,Qⱼ,Qᵢⱼ) = g(Pⱼ,Pᵢ,Qⱼ,Qᵢ,Qᵢⱼ)
g(Pᵢ,Pᵢ,Qᵢ,Qᵢ,Qᵢ) = Pᵢ
```
it is a more general form of `kij_mix`, where `kij_mix(f,P,Q) == pair_mix(g,P,Q)` is correct if:
```
f(Pᵢ,Pⱼ,Qᵢⱼ) = g(Pᵢ,Pⱼ,_,_,Qᵢⱼ)
```
"""
function pair_mix(f::F,P::ClapeyronParam,Q::ClapeyronParam) where F
    out = PairParam(P) #copy the input
    if Q isa SingleParameter
        vals = Diagonal(Q.vals)
        missingvals = out.ismissingvalues
        Q = PairParam(Q.name,Q.components,vals,true,out.ismissingvalues,Q.sources,Q.sourcecsvs)
    end

    return pair_mix!(f,out,Q)
end

function pair_mix!(f::F,out::PairParameter,Q::PairParameter) where F
    out_missing = out.ismissingvalues
    pair_mix!(f,out.values,Q.values,out_missing)
    #consider the two here:
    out_missing .= out_missing .& Q.ismissingvalues
    #but diagonals are all non-missing, by default:
     for i in diagind(out_missing)
        out_missing[i] = false
    end
    return out
end

function pair_mix!(f,p,q,B)
    N = LinearAlgebra.checksquare(p)
    for j ∈ 1:N
        p_j = p[j,j]
        q_j = q[j,j]
        for i ∈ 1:N
            if B[j,i]
                p_i = p[i,i]
                q_i = q[i,i]
                q_ji =  q[j,i]
                p_ji = f(p_j,p_i,q_j,q_i,q_ji)
                p[j,i] = p_ji
            end
        end
    end
    return p
end

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
function group_sum(groups::GroupParam, ::Nothing)
    v = groups.n_groups_cache
    return [sum(vi) for vi in v]
end

group_sum(groups::GroupParam) = group_sum(groups, nothing)

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
    _0 = zero(eltype(p))
    res = zeros(eltype(p),length(lcomps))
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
export kij_mix, pair_mix
export group_sum,group_pairmean
