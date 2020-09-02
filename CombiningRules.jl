# Lorentz-Berthelot rules

function combining_σ(σ)
    return sum([σ[i] for i in keys(σ)])/length(σ)
end

function combining_ϵ(ϵ)
    return prod([ϵ[i] for i in keys(ϵ)])^(1/length(ϵ))
end
