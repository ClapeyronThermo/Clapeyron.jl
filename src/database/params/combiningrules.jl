"""
    sigma_LorentzBerthelot(σ::ClapeyronParam,ζ::PairParam)::PairParam
    sigma_LorentzBerthelot(σ::ClapeyronParam)::PairParam

Combining rule for a single or pair parameter. returns a pair parameter with non diagonal entries equal to:
```
σᵢⱼ = (1 - ζᵢⱼ)*(σᵢ + σⱼ)/2
```
If `ζᵢⱼ` is not defined, the definition is reduced to a simple arithmetic mean:
```
σᵢⱼ = (σᵢ + σⱼ)/2
```
Ignores non-diagonal entries already set.

If a Single Parameter is passed as input, it will be converted to a Pair Parameter with `σᵢᵢ = σᵢ`.
"""
function sigma_LorentzBerthelot end

function sigma_LorentzBerthelot(sigma::SingleOrPair,zeta::PairParameter) 
    return sigma_LorentzBerthelot!(PairParam(sigma),zeta)
end

function sigma_LorentzBerthelot(sigma::SingleOrPair) 
    return sigma_LorentzBerthelot!(PairParam(sigma))
end

"""
    epsilon_LorentzBerthelot(ϵ::ClapeyronParam,k::PairParam)::PairParam
    epsilon_LorentzBerthelot(ϵ::ClapeyronParam)::PairParam

Combining rule for a single or pair parameter. returns a pair parameter with non diagonal entries equal to:
```
ϵᵢⱼ = (1 - kᵢⱼ)*√(ϵᵢϵⱼ)
```
If `kᵢⱼ` is not defined, the definition is reduced to a simple geometric mean:
```
ϵᵢⱼ = √(ϵᵢϵⱼ)
```
Ignores non-diagonal entries already set.

If a Single Parameter is passed as input, it will be converted to a Pair Parameter with `ϵᵢᵢ = ϵᵢ`.
"""
function epsilon_LorentzBerthelot end

function epsilon_LorentzBerthelot(epsilon::SingleOrPair, k::PairParameter)
    return epsilon_LorentzBerthelot!(PairParam(epsilon),k)
end

function epsilon_LorentzBerthelot(epsilon::SingleOrPair)
    return epsilon_LorentzBerthelot!(PairParam(epsilon))
end

"""
    epsilon_HudsenMcCoubrey(ϵ::ClapeyronParam,σ::PairParam)::PairParam
    epsilon_HudsenMcCoubrey(ϵ::ClapeyronParam)::PairParam

Combining rule for a single or pair parameter. returns a pair parameter with non diagonal entries equal to:
```
ϵᵢⱼ = √(ϵᵢϵⱼ)*(σᵢᵢ^3 * σⱼⱼ^3)/σᵢⱼ^6 
```
If `σᵢⱼ` is not defined, the definition is reduced to a simple geometric mean:
```
ϵᵢⱼ = √(ϵᵢϵⱼ)
```
Ignores non-diagonal entries already set.

If a Single Parameter is passed as input, it will be converted to a Pair Parameter with `ϵᵢᵢ = ϵᵢ`.
"""
function epsilon_HudsenMcCoubrey(epsilon::SingleOrPair, sigma::PairParameter)
    return epsilon_HudsenMcCoubrey!(PairParam(epsilon),sigma)
end

epsilon_HudsenMcCoubrey(epsilon) = epsilon_LorentzBerthelot(epsilon)

"""
    lambda_LorentzBerthelot(λ::ClapeyronParam,k::Real = 3)::PairParam

Combining rule for a single or pair parameter. returns a pair parameter with non diagonal entries equal to:
```
λᵢⱼ = k + √((λᵢᵢ - k)(λⱼⱼ - k))
```
with `k = 0` the definition is reduced to a simple geometric mean:
```
ϵᵢⱼ = √(ϵᵢϵⱼ)
```
Ignores non-diagonal entries already set.

If a Single Parameter is passed as input, it will be converted to a Pair Parameter with `λᵢᵢ = λᵢ`.
"""
function lambda_LorentzBerthelot(lambda::ClapeyronParam,k = 3)
    return lambda_LorentzBerthelot!(PairParam(lambda),k)
end

"""
    lambda_squarewell(λ::Union{PairParameter,SingleParameter},σ::PairParam)::PairParam

Combining rule for a single or pair parameter. returns a pair parameter with non diagonal entries equal to:
```
λᵢⱼ = (σᵢᵢλᵢᵢ + σⱼⱼλⱼⱼ)/(σᵢᵢ + σⱼⱼ)
```
Ignores non-diagonal entries already set.

If a Single Parameter is passed as input, it will be converted to a Pair Parameter with `λᵢᵢ = λᵢ`.
"""
function lambda_squarewell end

function lambda_squarewell(lambda::Union{PairParameter,SingleParameter}, sigma::Union{PairParameter,SingleParameter})
    return lambda_squarewell!(PairParam(lambda),sigma)
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
export kij_mix, pair_mix
export sigma_LorentzBerthelot
export epsilon_LorentzBerthelot
export epsilon_HudsenMcCoubrey
export lambda_LorentzBerthelot
export lambda_squarewell