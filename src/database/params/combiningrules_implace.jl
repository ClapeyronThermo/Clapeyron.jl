mix_mean(p_i,p_j,k=0) = 0.5*(p_i+p_j)*(1-k)  
mix_geomean(p_i,p_j,k=0) = sqrt(p_i*p_j)*(1-k) 
mix_powmean(p_i,p_j,k=0,n=2) =(1-k)*(0.5*(p_i^n + p_j^n))^(1/n)

##special lambda with custom k
function mix_lambda(λ_i,λ_j,k)
    return k + sqrt((λ_i - k) * (λ_j - k))
end

#general pair mix functions
mix_HudsenMcCoubrey(ϵᵢ,ϵⱼ,σᵢ,σⱼ,σᵢⱼ) = √(ϵᵢ*ϵⱼ)*(σᵢ^3 * σⱼ^3)/σᵢⱼ^6 
mix_lambda_squarewell(λᵢ,λⱼ,σᵢ,σⱼ,σᵢⱼ) = (σᵢ*λᵢ + σⱼ*λⱼ)/(σᵢ + σⱼ)

"""
    sigma_LorentzBerthelot!(σ::PairParam,ζ::PairParam)::PairParam
    sigma_LorentzBerthelot!(σ::PairParam)::PairParam

Combining rule for a pair parameter. returns a pair parameter with non diagonal entries equal to:
```
σᵢⱼ = (1 - ζᵢⱼ)*(σᵢ + σⱼ)/2
```
If `ζᵢⱼ` is not defined, the definition is reduced to a simple arithmetic mean:
```
σᵢⱼ = (σᵢ + σⱼ)/2
```
Ignores non-diagonal entries already set. and modifies `σ` implace.

For the out of place version, check [`sigma_LorentzBerthelot`](@ref)
"""
function sigma_LorentzBerthelot! end

sigma_LorentzBerthelot!(sigma::PairParameter,zeta::PairParameter) = kij_mix!(mix_mean,sigma,zeta)
sigma_LorentzBerthelot!(sigma::PairParameter) = kij_mix!(mix_mean,sigma,nothing)

"""
    epsilon_LorentzBerthelot!(ϵ::ClapeyronParam,k::PairParam)::PairParam
    epsilon_LorentzBerthelot!(ϵ::ClapeyronParam)::PairParam

Combining rule for a pair parameter. returns a pair parameter with non diagonal entries equal to:
```
ϵᵢⱼ = (1 - kᵢⱼ)*√(ϵᵢϵⱼ)
```
If `kᵢⱼ` is not defined, the definition is reduced to a simple geometric mean:
```
ϵᵢⱼ = √(ϵᵢϵⱼ)
```
Ignores non-diagonal entries already set. and modifies `ϵ` implace.

For the out of place version, check [`epsilon_LorentzBerthelot`](@ref)
"""
function epsilon_LorentzBerthelot! end

epsilon_LorentzBerthelot!(epsilon::ClapeyronParam, k::PairParameter) = kij_mix!(mix_geomean,epsilon,k)
epsilon_LorentzBerthelot!(epsilon::ClapeyronParam) = kij_mix!(mix_geomean,epsilon)

"""
    epsilon_HudsenMcCoubrey!(ϵ::PairParam,σ::PairParam)::PairParam
    epsilon_HudsenMcCoubrey!(ϵ::PairParam)::PairParam

Combining rule for a single or pair parameter. returns a pair parameter with non diagonal entries equal to:
```
ϵᵢⱼ = √(ϵᵢϵⱼ)*(σᵢᵢ^3 * σⱼⱼ^3)/σᵢⱼ^6 
```
If `σᵢⱼ` is not defined, the definition is reduced to a simple geometric mean:
```
ϵᵢⱼ = √(ϵᵢϵⱼ)
```
Ignores non-diagonal entries already set. and modifies `ϵ` implace.

For the out of place version, check [`epsilon_LorentzBerthelot`](@ref)
"""
function epsilon_HudsenMcCoubrey!(epsilon::PairParameter, sigma::PairParameter)
    return pair_mix!(mix_HudsenMcCoubrey,epsilon,sigma)
end

epsilon_HudsenMcCoubrey!(epsilon) = epsilon_LorentzBerthelot!(epsilon)

"""
    lambda_LorentzBerthelot!(λ::PairParam,k::Real = 3)::PairParam

Combining rule for a pair parameter. returns a pair parameter with non diagonal entries equal to:
```
λᵢⱼ = k + √((λᵢᵢ - k)(λⱼⱼ - k))
```
with `k = 0` the definition is reduced to a simple geometric mean:
```
ϵᵢⱼ = √(ϵᵢϵⱼ)
```
Ignores non-diagonal entries already set.

Ignores non-diagonal entries already set. and modifies `λ` implace.

For the out of place version, check [`lambda_LorentzBerthelot`](@ref)
"""
function lambda_LorentzBerthelot!(lambda::PairParam,k = 3)
    f(λi,λj,m) = mix_lambda(λi,λj,k) 
    return kij_mix!(f,lambda)
end

"""
    lambda_squarewell!(λ::PairParam,σ::PairParam)::PairParam

Combining rule for a pair parameter. returns a pair parameter with non diagonal entries equal to:
```
λᵢⱼ = (σᵢᵢλᵢᵢ + σⱼⱼλⱼⱼ)/(σᵢᵢ + σⱼⱼ)
```
Ignores non-diagonal entries already set.

Ignores non-diagonal entries already set. and modifies `λ` implace.

For the out of place version, check [`lambda_squarewell`](@ref)
"""
function lambda_squarewell!(lambda::PairParameter, sigma::PairParameter)
    return pair_mix!(mix_lambda_squarewell,lambda,sigma)
end

export sigma_LorentzBerthelot!
export epsilon_LorentzBerthelot!
export epsilon_HudsenMcCoubrey!
export lambda_LorentzBerthelot!
export lambda_squarewell!