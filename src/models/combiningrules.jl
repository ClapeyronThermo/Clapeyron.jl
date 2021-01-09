function combining_sigma(sigma::OpenSAFTParam; rules="Lorentz-Berthelot")
    sigma = PairParam(sigma)
    σ = sigma.values
    for i ∈ 1:length(sigma.components), j ∈ 1:length(sigma.components)
        if sigma.ismissingvalues[i,j]
            if rules == "Lorentz-Berthelot"
                σ[i,j] = (σ[i,i] + σ[j,j]) / 2
            end
        end
    end
    return sigma
end

function combining_epsilon(epsilon::OpenSAFTParam, combiningparam::PairParam; rules="Lorentz-Berthelot")
    epsilon = PairParam(epsilon)
    ϵ = epsilon.values
    if rules == "Lorentz-Berthelot"
        k = combiningparam
        k_ = k.values
    elseif rules == "Hudsen-McCoubrey"
        sigma = combiningparam
        σ = sigma.values
    end
    for i ∈ 1:length(epsilon.components), j ∈ 1:length(epsilon.components)
        if epsilon.ismissingvalues[i,j]
            if rules == "Lorentz-Berthelot"
                ϵ[i,j] = (1 - k_[i,j]) * sqrt(ϵ[i,i] * ϵ[j,j])
            elseif rules == "Hudsen-McCoubrey"
                ϵ[i,j] = (σ[i,i]^3 * σ[j,j]^3) / σ[i,j]^6 * sqrt(ϵ[i,i] * ϵ[j,j])
            end
        end
    end
    return epsilon
end

function combining_lambda(lambda::OpenSAFTParam; rules="Lorentz-Berthelot")
    lambda = PairParam(lambda)
    λ = lambda.values
    for i ∈ 1:length(lambda.components), j ∈ 1:length(lambda.components)
        if lambda.ismissingvalues[i,j]
            if rules == "Lorentz-Berthelot"
                λ[i,j] = 3 + sqrt((λ[i,i] - 3) * (λ[j,j] - 3))
            end
        end
    end
    return lambda
end

function combining_lambda(lambda::OpenSAFTParam, combiningparam::PairParam; rules="Square-Well")
    lambda = PairParam(lambda)
    sigma = combiningparam
    λ = lambda.values
    σ = sigma.values
    for i ∈ 1:length(lambda.components), j ∈ 1:length(lambda.components)
        if lambda.ismissingvalues[i,j]
            if rules == "Square-Well"
                λ[i,j] = (λ[i,i]*σ[i,i] + λ[j,j]*σ[j,j]) / (σ[i,i] + σ[j,j])
            end
        end
    end
    return lambda
end
