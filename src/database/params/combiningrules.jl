### base functions to use:
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

export sigma_LorentzBerthelot
export epsilon_LorentzBerthelot
export epsilon_HudsenMcCoubrey
export lambda_LorentzBerthelot
export lambda_squarewell