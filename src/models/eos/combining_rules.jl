function combining_sigma(sigma::Dict; rules="Lorentz-Berthelot")
    components = [collect(i)[1] for i in keys(sigma)]
    pairs = vcat([Set([i]) for i in components], [Set(i) for i in collect(combinations(components, 2))])
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

function combining_epsilon(epsilon::Dict, sigma::Dict, k::Dict; rules_k="Lorentz-Berthelot-mod", rules_no_k="Lorentz-Berthelot")
    components = [collect(i)[1] for i in keys(epsilon)]
    pairs = vcat([Set([i]) for i in components], [Set(i) for i in collect(combinations(components, 2))])
    combined_epsilons = Dict{Set{String},Float64}()
    for pair in pairs
        if !haskey(epsilon, pair)
            if haskey(k, pair)
                if rules_k == "Lorentz-Berthelot-mod"
                    combined_epsilon = (1-k[pair])*sqrt(prod(epsilon[Set([i])] for i in pair))
                end
            else
                if rules_no_k == "Lorentz-Berthelot"
                    combined_epsilon = sqrt(prod(epsilon[Set([i])] for i in pair))
                elseif rules_no_k == "Hudson-McCoubrey"
                    combined_epsilon = 2^3*sqrt(prod(sigma[Set([i])]^3 for i in pair))/sum(sigma[Set([i])] for i in pair)^3*sqrt(prod(epsilon[Set([i])] for i in pair))
                end
            end
            push!(combined_epsilons, pair => combined_epsilon)
        end
    end
    return combined_epsilons
end

function combining_lambda_Mie(lambda::Dict; rules="Lorentz-Berthelot")
    components = [collect(i)[1] for i in keys(lambda)]
    pairs = vcat([Set([i]) for i in components], [Set(i) for i in collect(combinations(components, 2))])
    combined_lambdas = Dict{Set{String},Float64}()
    for pair in pairs
        if !haskey(lambda, pair)
            if rules == "Lorentz-Berthelot"
                combined_lambda = 3+sqrt(prod(lambda[Set([i])]-3 for i in pair))
            end
            push!(combined_lambdas, pair => combined_lambda)
        end
    end
    return combined_lambdas
end

function combining_lambda_SW(lambda::Dict,sigma::Dict)
    components = [collect(i)[1] for i in keys(lambda)]
    pairs = vcat([Set([i]) for i in components], [Set(i) for i in collect(combinations(components, 2))])
    combined_lambdas = Dict{Set{String},Float64}()
    for pair in pairs
        if !haskey(lambda, pair)
            combined_lambda = sum(lambda[Set([i])]*sigma[Set([i])] for i in pair)/sum(sigma[Set([i])] for i in pair)
            push!(combined_lambdas, pair => combined_lambda)
        end
    end
    return combined_lambdas
end
