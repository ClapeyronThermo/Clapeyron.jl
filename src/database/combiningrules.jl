include("combiningrules/base.jl")
include("combiningrules/implace.jl")
include("combiningrules/outplace.jl")
include("combiningrules/group.jl")
include("combiningrules/assoc.jl")
"""
    sigma_LorentzBerthelot(σ::SingleOrPair,ζ::PairParam)::PairParam
    sigma_LorentzBerthelot(σ::SingleOrPair)::PairParam

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

sigma_LorentzBerthelot(sigma::SingleOrPair,zeta = nothing) = kij_mix(mix_mean,sigma,zeta)

sigma_LorentzBerthelot!(sigma::PairParameter,zeta = nothing) = kij_mix!(mix_mean,sigma,zeta)


"""
    epsilon_LorentzBerthelot(ϵ::SingleOrPair,k::PairParam)::PairParam
    epsilon_LorentzBerthelot(ϵ::SingleOrPair)::PairParam

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

epsilon_LorentzBerthelot(epsilon::SingleOrPair, k = nothing) = kij_mix(mix_geomean,epsilon,k)

epsilon_LorentzBerthelot!(epsilon::PairParameter, k = nothing) = kij_mix!(mix_geomean,epsilon,k)


"""
    epsilon_HudsenMcCoubrey(ϵ::SingleOrPair,σ::PairParam)::PairParam
    epsilon_HudsenMcCoubrey(ϵ::SingleOrPair)::PairParam

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
function epsilon_HudsenMcCoubrey(epsilon::SingleOrPair, sigma::PairParameter;k = nothing)
    return pair_mix(mix_HudsenMcCoubrey,epsilon,sigma)
end

epsilon_HudsenMcCoubrey(epsilon) = epsilon_LorentzBerthelot(epsilon)
epsilon_HudsenMcCoubrey(epsilon,::Nothing) = epsilon_LorentzBerthelot(epsilon)

function epsilon_HudsenMcCoubrey!(epsilon::PairParameter, sigma::PairParameter)
    return pair_mix!(mix_HudsenMcCoubrey,epsilon,sigma)
end

epsilon_HudsenMcCoubrey(epsilon::PairParameter) = epsilon_LorentzBerthelot!(epsilon)


function lambda_LorentzBerthelot!(lambda::PairParameter,k = 3)
    f(λi,λj,m) = mix_lambda(λi,λj,k) 
    return kij_mix!(f,lambda)
end

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
function lambda_LorentzBerthelot(lambda::SingleOrPair,k = 3)
    return return lambda_LorentzBerthelot!(PairParam(lambda),k)
end


"""
    lambda_squarewell(λ::SingleOrPair,σ::PairParam)::PairParam

Combining rule for a single or pair parameter. returns a pair parameter with non diagonal entries equal to:
```
λᵢⱼ = (σᵢᵢλᵢᵢ + σⱼⱼλⱼⱼ)/(σᵢᵢ + σⱼⱼ)
```
Ignores non-diagonal entries already set.

If a Single Parameter is passed as input, it will be converted to a Pair Parameter with `λᵢᵢ = λᵢ`.
"""
function lambda_squarewell end

function lambda_squarewell(lambda::SingleOrPair, sigma::Union{PairParameter,SingleParameter})
    return pair_mix(mix_lambda_squarewell,lambda,sigma)
end

function lambda_squarewell!(lambda::PairParameter, sigma::Union{PairParameter,SingleParameter})
    return pair_mix!(mix_lambda_squarewell,lambda,sigma)
end

export kij_mix, pair_mix
export sigma_LorentzBerthelot
export epsilon_LorentzBerthelot
export epsilon_HudsenMcCoubrey
export lambda_LorentzBerthelot
export lambda_squarewell
export group_sum,group_pairmean