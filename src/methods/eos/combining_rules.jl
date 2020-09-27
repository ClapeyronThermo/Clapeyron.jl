function combining_sigma(sigma::Dict; rules="Lorentz-Berthelot")
    components = [collect(i)[1] for i in keys(sigma)]
    pairs = [Set(i) for i in collect(combinations(components, 2))]
    combined_sigmas = Dict{Set{String},Float64}()
    for pair in pairs
        if !haskey(sigma, pair)
            if rules == "Lorentz-Berthelot"
                combined_sigma = (sum(sigma[Set([i])] for i in pair))/2
            end
            push!(combined_sigmas, pair => combined_sigma)
        end
    end
    return combined_sigmas
end

function combining_epsilon(epsilon::Dict, k::Dict; rules_k="Hudsen-McCoubrey", rules_no_k="Lorentz-Berthelot")
    components = [collect(i)[1] for i in keys(epsilon)]
    pairs = [Set(i) for i in collect(combinations(components, 2))]
    combined_epsilons = Dict{Set{String},Float64}()
    for pair in pairs
        if !haskey(epsilon, pair)
            if haskey(k, pair)
                if rules_k == "Hudsen-McCoubrey"
                    combined_epsilon = (1-k[pair])*sqrt(prod(epsilon[Set([i])] for i in pair))
                end
            else
                if rules_no_k == "Lorentz-Berthelot"
                    combined_epsilon = sqrt(prod(epsilon[Set([i])] for i in pair))
                end
            end
            push!(combined_epsilons, pair => combined_epsilon) 
        end
    end
    return combined_epsilons
end
