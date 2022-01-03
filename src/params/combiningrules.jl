export sigma_LorentzBerthelot_mod
function sigma_LorentzBerthelot_mod(sigma::ClapeyronParam,zeta::PairParameter)
    sigma = PairParam(sigma)
    σ = sigma.values
    zeta_ = zeta.values
    for i ∈ 1:length(sigma.components), j ∈ 1:length(sigma.components)
        if sigma.ismissingvalues[i,j]
            σ[i,j] = (1-zeta_[i,j])*(σ[i,i] + σ[j,j]) / 2
        end
    end
    return sigma
end

export sigma_LorentzBerthelot
function sigma_LorentzBerthelot(sigma::ClapeyronParam)
    sigma = PairParam(sigma)
    σ = sigma.values
    for i ∈ 1:length(sigma.components), j ∈ 1:length(sigma.components)
        if sigma.ismissingvalues[i,j]
            σ[i,j] = (σ[i,i] + σ[j,j]) / 2
        end
    end
    sigma.ismissingvalues .= false
    return sigma
end

export epsilon_LorentzBerthelot
function epsilon_LorentzBerthelot(epsilon::ClapeyronParam, k::PairParameter)
    epsilon = PairParam(epsilon)
    ϵ = epsilon.values
    k_ = k.values
    for i ∈ 1:length(epsilon.components), j ∈ 1:length(epsilon.components)
        if epsilon.ismissingvalues[i,j]
            ϵ[i,j] = (1 - k_[i,j]) * sqrt(ϵ[i,i] * ϵ[j,j])
        end
    end
    epsilon.ismissingvalues .= false
    return epsilon
end

export epsilon_HudsenMcCoubrey
function epsilon_HudsenMcCoubrey(epsilon::ClapeyronParam, sigma::PairParameter)
    epsilon = PairParam(epsilon)
    ϵ = epsilon.values
    σ = sigma.values
    for i ∈ 1:length(epsilon.components), j ∈ 1:length(epsilon.components)
        if epsilon.ismissingvalues[i,j]
            ϵ[i,j] = (σ[i,i]^3 * σ[j,j]^3) / σ[i,j]^6 * sqrt(ϵ[i,i] * ϵ[j,j])
        end
    end
    epsilon.ismissingvalues .= false
    return epsilon
end

export lambda_LorentzBerthelot
function lambda_LorentzBerthelot(lambda::ClapeyronParam)
    lambda = PairParam(lambda)
    λ = lambda.values
    for i ∈ 1:length(lambda.components), j ∈ 1:length(lambda.components)
        if lambda.ismissingvalues[i,j]
            λ[i,j] = 3 + sqrt((λ[i,i] - 3) * (λ[j,j] - 3))
        end
    end
    lambda.ismissingvalues .= false
    return lambda
end

export lambda_squarewell
function lambda_squarewell(lambda::ClapeyronParam, sigma::PairParameter)
    lambda = PairParam(lambda)
    λ = lambda.values
    σ = sigma.values
    for i ∈ 1:length(lambda.components), j ∈ 1:length(lambda.components)
        if lambda.ismissingvalues[i,j]
            λ[i,j] = (λ[i,i]*σ[i,i] + λ[j,j]*σ[j,j]) / (σ[i,i] + σ[j,j])
        end
    end
    lambda.ismissingvalues .= false
    return lambda
end

"""
    mixing_rule_quad(op, x, p)


returns an efficient implementation of:
sum(x[i]*x[j]*op(p[i],p[j]) for i in @comps for j in @comps)`
where `op` is a function of two arguments that satisfies op(p[i],p[i]) = p[i]

# example
```julia-repl
julia> prop = [200,600];
julia> x0 = [0.3,0.7];
julia> propmix1 = Clapeyron.mixing_rule_quad((x,y)->0.5*(x+y),x0,prop);
julia> propmix2 = sum(x0[i]*x0[j]*(prop[i]+prop[j])*0.5 for i in 1:2 for j in 1:2);
julia> propmix1 ≈ propmix2;
true
```
"""
function mixing_rule_quad(op, x, p)
    N = length(x)
    @boundscheck checkbounds(x, N)
    @inbounds begin
        res1 = zero(eltype(x))
        for i = 1:N
            res1 += p[i] * x[i]^2
            for j = 1:i - 1
                res1 += 2 * x[i] * x[j] * op(p[i], p[j])
            end
        end
    end
    return res1
end
