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

export sigma_LorentzBerthelot
export epsilon_LorentzBerthelot
export epsilon_HudsenMcCoubrey
export lambda_LorentzBerthelot
export lambda_squarewell