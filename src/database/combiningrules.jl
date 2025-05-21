include("combiningrules/base.jl")
include("combiningrules/implace.jl")
include("combiningrules/outplace.jl")
include("combiningrules/group.jl")
include("combiningrules/assoc.jl")
include("combiningrules/get_k.jl")

"""
    sigma_LorentzBerthelot(σ::SingleOrPair,ζ)::PairParam
    sigma_LorentzBerthelot(σ::SingleOrPair)::PairParam
    sigma_LorentzBerthelot(σ::Union{AbstractVector,AbstractMatrix},ζ)::AbstractMatrix
    sigma_LorentzBerthelot(σ::Union{AbstractVector,AbstractMatrix})::AbstractMatrix

Combining rule for a single or pair parameter.
Returns a pair parameter with non diagonal entries equal to:

```
σᵢⱼ = (1 - ζᵢⱼ)*(σᵢ + σⱼ)/2
```

If `ζᵢⱼ` is not defined, the definition is reduced to a simple arithmetic mean:

```
σᵢⱼ = (σᵢ + σⱼ)/2
```

Ignores non-diagonal entries already set.

If a Single Parameter (or vector) is passed as input, it will be converted to a Pair Parameter with `σᵢᵢ = σᵢ`.
"""
sigma_LorentzBerthelot(sigma,zeta = nothing) = kij_mix(mix_mean,sigma,zeta)

"""
    sigma_LorentzBerthelot!(σ::PairParameter,ζ)::PairParam
    sigma_LorentzBerthelot!(σ::AbstractMatrix,ζ)::AbstractMatrix

Combining rule for a single or pair parameter.
Returns a pair parameter with non diagonal entries equal to:

```
σᵢⱼ = (1 - ζᵢⱼ)*(σᵢ + σⱼ)/2
```

If `ζᵢⱼ` is not defined, the definition is reduced to a simple arithmetic mean:

```
σᵢⱼ = (σᵢ + σⱼ)/2
```

The method overwrites the entries in `σ`, with the exception of diagonal entries.
"""
sigma_LorentzBerthelot!(sigma,zeta = nothing) = kij_mix!(mix_mean,sigma,zeta)

"""
    epsilon_LorentzBerthelot(ϵ::SingleOrPair,k)::PairParam
    epsilon_LorentzBerthelot(ϵ::SingleOrPair)::PairParam
    epsilon_LorentzBerthelot(ϵ::Union{AbstractVector,AbstractMatrix},k)::AbstractMatrix
    epsilon_LorentzBerthelot(ϵ::Union{AbstractVector,AbstractMatrix})::AbstractMatrix

Combining rule for a single or pair parameter.
Returns a pair parameter with non diagonal entries equal to:

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
epsilon_LorentzBerthelot(epsilon, k = nothing) = kij_mix(mix_geomean,epsilon,k)

"""
    epsilon_LorentzBerthelot!(ϵ::PairParameter,k)::PairParam
    epsilon_LorentzBerthelot!(ϵ::PairParameter)::PairParam
    epsilon_LorentzBerthelot!(ϵ::AbstractMatrix,k)::AbstractMatrix
    epsilon_LorentzBerthelot!(ϵ::AbstractMatrix)::AbstractMatrix

Combining rule for a single or pair parameter.
Returns a pair parameter with non diagonal entries equal to:

```
ϵᵢⱼ = (1 - kᵢⱼ)*√(ϵᵢϵⱼ)
```

If `kᵢⱼ` is not defined, the definition is reduced to a simple geometric mean:

```
ϵᵢⱼ = √(ϵᵢϵⱼ)
```

The method overwrites the entries in `ϵ`, with the exception of diagonal entries.
"""
epsilon_LorentzBerthelot!(epsilon, k = nothing) = kij_mix!(mix_geomean,epsilon,k)

"""
    epsilon_HudsenMcCoubrey(ϵ::SingleOrPair,σ)::PairParam
    epsilon_HudsenMcCoubrey(ϵ::SingleOrPair)::PairParam
    epsilon_HudsenMcCoubrey(ϵ::Union{AbstractVector,AbstractMatrix},σ)::AbstractMatrix
    epsilon_HudsenMcCoubrey(ϵ::Union{AbstractVector,AbstractMatrix})::AbstractMatrix

Combining rule for a single or pair parameter.
Returns a pair parameter with non diagonal entries equal to:

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
function epsilon_HudsenMcCoubrey end
epsilon_HudsenMcCoubrey(epsilon, sigma) = pair_mix(mix_HudsenMcCoubrey,epsilon,sigma)
epsilon_HudsenMcCoubrey(epsilon) = epsilon_LorentzBerthelot(epsilon)
epsilon_HudsenMcCoubrey(epsilon,::Nothing) = epsilon_HudsenMcCoubrey(epsilon)

"""
    epsilon_HudsenMcCoubrey!(ϵ::PairParameter,σ)::PairParam
    epsilon_HudsenMcCoubrey!(ϵ::PairParameter)::PairParam
    epsilon_HudsenMcCoubrey!(ϵ::AbstractMatrix,σ)::AbstractMatrix
    epsilon_HudsenMcCoubrey!(ϵ::AbstractMatrix)::AbstractMatrix

Combining rule for a single or pair parameter.
Returns a pair parameter with non diagonal entries equal to:

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
function epsilon_HudsenMcCoubrey! end
epsilon_HudsenMcCoubrey!(epsilon, sigma) = pair_mix!(mix_HudsenMcCoubrey,epsilon,sigma)
epsilon_HudsenMcCoubrey!(epsilon) = epsilon_LorentzBerthelot!(epsilon)
epsilon_HudsenMcCoubrey!(epsilon,::Nothing) = epsilon_HudsenMcCoubrey!(epsilon)

"""
    epsilon_HudsenMcCoubreysqrt(ϵ::SingleOrPair,σ)::PairParam
    epsilon_HudsenMcCoubreysqrt(ϵ::SingleOrPair)::PairParam
    epsilon_HudsenMcCoubreysqrt(ϵ::Union{AbstractVector,AbstractMatrix},σ)::AbstractMatrix
    epsilon_HudsenMcCoubreysqrt(ϵ::Union{AbstractVector,AbstractMatrix})::AbstractMatrix

Combining rule for a single or pair parameter. returns a pair parameter with non diagonal entries equal to:
```
ϵᵢⱼ = √(ϵᵢϵⱼ * σᵢᵢ^3 * σⱼⱼ^3)/σᵢⱼ^3
```
If `σᵢⱼ` is not defined, the definition is reduced to a simple geometric mean:
```
ϵᵢⱼ = √(ϵᵢϵⱼ)
```
Ignores non-diagonal entries already set.

If a Single Parameter is passed as input, it will be converted to a Pair Parameter with `ϵᵢᵢ = ϵᵢ`.
"""
function epsilon_HudsenMcCoubreysqrt end
epsilon_HudsenMcCoubreysqrt(epsilon, sigma) = pair_mix(mix_HudsenMcCoubreysqrt,epsilon,sigma)
epsilon_HudsenMcCoubreysqrt(epsilon) = epsilon_LorentzBerthelot(epsilon)
epsilon_HudsenMcCoubreysqrt(epsilon,::Nothing) = epsilon_LorentzBerthelot(epsilon)

function epsilon_HudsenMcCoubreysqrt! end
epsilon_HudsenMcCoubreysqrt!(epsilon) = epsilon_LorentzBerthelot!(epsilon)
epsilon_HudsenMcCoubreysqrt!(epsilon, sigma) = pair_mix!(mix_HudsenMcCoubreysqrt,epsilon,sigma)
epsilon_HudsenMcCoubreysqrt!(epsilon, ::Nothing) = epsilon_HudsenMcCoubreysqrt!(epsilon)

"""
    lambda_LorentzBerthelot(λ::SingleOrPair,k = 3)::PairParam
    lambda_LorentzBerthelot(λ::Union{AbstractVector,AbstractMatrix},k = 3)::AbstractMatrix

Combining rule for a single or pair parameter. returns a pair parameter with non diagonal entries equal to:
```
λᵢⱼ = k + √((λᵢᵢ - k)(λⱼⱼ - k))
```
with `k = 0` the definition is reduced to a simple geometric mean:
```
λᵢⱼ = √(λᵢλⱼ)
```
Ignores non-diagonal entries already set.

If a Single Parameter is passed as input, it will be converted to a Pair Parameter with `λᵢᵢ = λᵢ`.
"""
lambda_LorentzBerthelot(lambda,k = 3) = kij_mix(MixLambda(k),lambda)
lambda_LorentzBerthelot!(lambda::PairParameter,k = 3) = kij_mix!(MixLambda(k),lambda)

"""
    lambda_squarewell(λ::SingleOrPair,σ)::PairParam
    lambda_squarewell(λ::Union{AbstractVector,AbstractMatrix},σ)::AbstractMatrix

Combining rule for a single or pair parameter. returns a pair parameter with non diagonal entries equal to:
```
λᵢⱼ = (σᵢᵢλᵢᵢ + σⱼⱼλⱼⱼ)/(σᵢᵢ + σⱼⱼ)
```
Ignores non-diagonal entries already set.

If a Single Parameter is passed as input, it will be converted to a Pair Parameter with `λᵢᵢ = λᵢ`.
"""
lambda_squarewell(lambda, sigma) = pair_mix(mix_lambda_squarewell,lambda,sigma)
lambda_squarewell!(lambda, sigma) = pair_mix!(mix_lambda_squarewell,lambda,sigma)

function mirror_pair!(param::PairParameter,f = identity)
    mirror_pair!(param.values,param.ismissingvalues,f)
    return param
end

"""
    mirror_pair(param::PairParam,f = identity)

performs an operation `f` over the indices of `p` such as `p[j,i] = f(p[i,j])`. by default, `f = identity` (a symmetric matrix). 
One key difference is that it sets the `ismissingvalues` field for each modified index to `false`
"""
mirror_pair(param::PairParameter,f = identity) = mirror_pair!(deepcopy(param),f)

export kij_mix, kij_mix!
export pair_mix, pair_mix!
export sigma_LorentzBerthelot, sigma_LorentzBerthelot!
export epsilon_LorentzBerthelot, epsilon_LorentzBerthelot!
export epsilon_HudsenMcCoubrey, epsilon_HudsenMcCoubrey!
export epsilon_HudsenMcCoubreysqrt, epsilon_HudsenMcCoubreysqrt!
export lambda_LorentzBerthelot, lambda_LorentzBerthelot!
export lambda_squarewell, lambda_squarewell!
export group_sum, group_sum!
export group_pairmean, group_pairmean!
