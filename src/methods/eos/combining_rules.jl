function combining_σ(σ, rules="Lorentz-Berthelot")
    if rules = "Lorentz-Berthelot"
        return sum([σ[i] for i in keys(σ)])/length(σ)
    end
end

function combining_ϵ(ϵ, rules="Lorentz-Berthelot")
    if rules = "Lorentz-Berthelot"
        return prod([ϵ[i] for i in keys(ϵ)])^(1/length(ϵ))
    end
end
