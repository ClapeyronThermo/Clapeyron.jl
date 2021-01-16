function sigma_LorentzBerthelot(sigma::OpenSAFTParam)
    sigma = PairParam(sigma)
    σ = sigma.values
    for i ∈ 1:length(sigma.components), j ∈ 1:length(sigma.components)
        if sigma.ismissingvalues[i,j]
            σ[i,j] = (σ[i,i] + σ[j,j]) / 2
        end
    end
    return sigma
end

function epsilon_LorentzBerthelot(epsilon::OpenSAFTParam, k::PairParam)
    epsilon = PairParam(epsilon)
    ϵ = epsilon.values
    k_ = k.values
    for i ∈ 1:length(epsilon.components), j ∈ 1:length(epsilon.components)
        if epsilon.ismissingvalues[i,j]
            ϵ[i,j] = (1 - k_[i,j]) * sqrt(ϵ[i,i] * ϵ[j,j])
        end
    end
    return epsilon
end

function epsilon_HudsenMcCoubrey(epsilon::OpenSAFTParam, sigma::PairParam)
    epsilon = PairParam(epsilon)
    ϵ = epsilon.values
    σ = sigma.values
    for i ∈ 1:length(epsilon.components), j ∈ 1:length(epsilon.components)
        if epsilon.ismissingvalues[i,j]
            ϵ[i,j] = (σ[i,i]^3 * σ[j,j]^3) / σ[i,j]^6 * sqrt(ϵ[i,i] * ϵ[j,j])
        end
    end
    return epsilon
end

function lambda_LorentzBerthelot(lambda::OpenSAFTParam)
    lambda = PairParam(lambda)
    λ = lambda.values
    for i ∈ 1:length(lambda.components), j ∈ 1:length(lambda.components)
        if lambda.ismissingvalues[i,j]
            λ[i,j] = 3 + sqrt((λ[i,i] - 3) * (λ[j,j] - 3))
        end
    end
    return lambda
end

function lambda_squarewell(lambda::OpenSAFTParam, sigma::PairParam)
    lambda = PairParam(lambda)
    λ = lambda.values
    σ = sigma.values
    for i ∈ 1:length(lambda.components), j ∈ 1:length(lambda.components)
        if lambda.ismissingvalues[i,j]
            λ[i,j] = (λ[i,i]*σ[i,i] + λ[j,j]*σ[j,j]) / (σ[i,i] + σ[j,j])
        end
    end
    return lambda
end
