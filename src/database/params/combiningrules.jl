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
function kij_mix(f::F,param::ClapeyronParam,K = nothing) where F
    N = length(param.components)
    
    out = PairParam(param)
    p = out.values
    if K === nothing
        k = FillArrays.Zeros(N,N)
        missingK = FillArrays.Fill(true,N,N)
    else
        k = K.values
        missingK = K.ismissingvalues
    end

    kij_mix!(f,p,k,out.ismissingvalues)
    #should consider the two.
    out.ismissingvalues .= out.ismissingvalues .& missingK
        
    #but diagonals are all non-missing, by default:
    for i in 1:N
        out.ismissingvalues[i,i] = false
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
### functions to use:
mix_mean(p_i,p_j,k=0) = 0.5*(p_i+p_j)*(1-k)  
mix_geomean(p_i,p_j,k=0) = sqrt(p_i*p_j)*(1-k) 
mix_powmean(p_i,p_j,k=0,n=2) =(1-k)*(0.5*(p_i^n + p_j^n))^(1/n)

##special lambda with custom k
function mix_lambda(λ_i,λ_j,k)
    return k + sqrt((λ_i - k) * (λ_j - k))
end
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

sigma_LorentzBerthelot(sigma::ClapeyronParam,zeta::PairParameter) = kij_mix(mix_mean,sigma,zeta)
sigma_LorentzBerthelot(sigma::ClapeyronParam) = kij_mix(mix_mean,sigma)

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

epsilon_LorentzBerthelot(epsilon::ClapeyronParam, k::PairParameter) = kij_mix(mix_geomean,epsilon,k)
epsilon_LorentzBerthelot(epsilon::ClapeyronParam) = kij_mix(mix_geomean,epsilon)

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
    out = PairParam(P)
    Q isa PairParameter || (Q = PairParam(Q))
    p = out.values
    q = Q.values
    missingP = out.ismissingvalues
    missingQ = Q.ismissingvalues
    
    pair_mix!(f,p,q,out.ismissingvalues)
    #consider the two here:
    out.ismissingvalues .= missingP .& missingQ
    #but diagonals are all non-missing, by default:
     for i in 1:size(out.ismissingvalues,1)
        out.ismissingvalues[i,i] = false
    end
    #out.ismissingvalues .= false
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

mix_HudsenMcCoubrey(ϵᵢ,ϵⱼ,σᵢ,σⱼ,σᵢⱼ) = √(ϵᵢ*ϵⱼ)*(σᵢ^3 * σⱼ^3)/σᵢⱼ^6 
mix_lambda_squarewell(λᵢ,λⱼ,σᵢ,σⱼ,σᵢⱼ) = (σᵢ*λᵢ + σⱼ*λⱼ)/(σᵢ + σⱼ)

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
function epsilon_HudsenMcCoubrey(epsilon::ClapeyronParam, sigma::PairParameter)
    return pair_mix(mix_HudsenMcCoubrey,epsilon,sigma)
end

epsilon_HudsenMcCoubrey(epsilon) = epsilon_LorentzBerthelot(epsilon)

"""
    lambda_LorentzBerthelot(λ::ClapeyronParam,k = 3)::PairParam

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
    f(λi,λj,m) = mix_lambda(λi,λj,k) 
    return kij_mix(f,lambda)
end

"""
    lambda_squarewell(λ::ClapeyronParam,σ::PairParam)::PairParam

Combining rule for a single or pair parameter. returns a pair parameter with non diagonal entries equal to:
```
λᵢⱼ = (σᵢᵢλᵢᵢ + σⱼⱼλⱼⱼ)/(σᵢᵢ + σⱼⱼ)
```
Ignores non-diagonal entries already set.

If a Single Parameter is passed as input, it will be converted to a Pair Parameter with `λᵢᵢ = λᵢ`.
"""
function lambda_squarewell end

function lambda_squarewell(lambda::ClapeyronParam, sigma::Union{PairParameter,SingleParameter})
    return pair_mix(mix_lambda_squarewell,lambda,sigma)
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
export sigma_LorentzBerthelot
export epsilon_LorentzBerthelot
export epsilon_HudsenMcCoubrey
export lambda_LorentzBerthelot
export lambda_squarewell
export group_sum,group_pairmean