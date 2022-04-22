abstract type ePCSAFTModel <: PCSAFTModel end
@newmodel ePCSAFT ePCSAFTModel PCSAFTParam

export ePCSAFT
function ePCSAFT(components; idealmodel=BasicIdeal, userlocations=String[], ideal_userlocations=String[], verbose=false)
    params,sites = getparams(components, ["SAFT/ePCSAFT"]; userlocations=userlocations, verbose=verbose)
    
    segment = params["m"]
    k = params["k"]
    Mw = params["Mw"]
    params["sigma"].values .*= 1E-10
    sigma = sigma_LorentzBerthelot(params["sigma"])
    epsilon = epsilon_LorentzBerthelot(params["epsilon"], k)
    epsilon_assoc = params["epsilon_assoc"]
    bondvol = params["bondvol"]
    packagedparams = PCSAFTParam(Mw, segment, sigma, epsilon, epsilon_assoc, bondvol)
    references = ["10.1021/ie020753p"]

    model = ePCSAFT(packagedparams, sites, idealmodel; ideal_userlocations=ideal_userlocations, references=references, verbose=verbose)
    return model
end

function d(model::ePCSAFTModel, V, T, z)
    return d.(model, V, T, z, @comps)
end

function d(model::ePCSAFTModel, V, T, z, i)
    ϵii = model.params.epsilon.diagvalues[i]

    if model.components[i]=="water"
        σii = (2.792700+10.1100*exp(-0.01775*T)-1.41700*exp(-0.01146*T))*1e-10
        return σii * (1 - 0.12exp(-3ϵii/T))
    elseif occursin("+",model.components[i]) || occursin("-",model.components[i])
        σii = model.params.sigma.diagvalues[i]
        return σii*(1-0.12)
    else
        σii = model.params.sigma.diagvalues[i]
        return σii * (1 - 0.12exp(-3ϵii/T))
    end
end

function m2ϵσ3(model::ePCSAFTModel, V, T, z)
    x = z/∑(z)
    m = model.params.segment.values
    σ = model.params.sigma.values
    ϵ = model.params.epsilon.values
    A = 0.
    B = 0.
    for i ∈ @comps
        for j ∈ @comps
            if model.components[i]=="water" && model.components[j]=="water"
                σij = (2.792700+10.1100*exp(-0.01775*T)-1.41700*exp(-0.01146*T))*1e-10
            elseif model.components[i]=="water" && model.components[j]!="water"
                σij = ((2.792700+10.1100*exp(-0.01775*T)-1.41700*exp(-0.01146*T))*1e-10+model.params.sigma.values[j,j])/2
            elseif model.components[j]=="water" && model.components[i]!="water"
                σij = ((2.792700+10.1100*exp(-0.01775*T)-1.41700*exp(-0.01146*T))*1e-10+model.params.sigma.values[i,i])/2
            else
                σij = model.params.sigma.values[i,j]
            end
            if occursin("+",model.components[i]) && occursin("-",model.components[j])
                A+=0.
                B+=0.
            elseif occursin("-",model.components[i]) && occursin("+",model.components[j])
                A+=0.
                B+=0.
            else
                A+=x[i]*x[j]*m[i]*m[j]* (ϵ[i,j]*(1)/T)^1 *σij^3
                B+=x[i]*x[j]*m[i]*m[j]* (ϵ[i,j]*(1)/T)^2 *σij^3
            end
        end
    end
    return A, B
end

# function Δ(model::ePCSAFTModel, V, T, z, i, j, a, b)
#     ϵ_assoc = model.params.epsilon_assoc.values
#     κ = model.params.bondvol.values
    
#     gij = @f(g_hs,i,j)
#     return gij*σij^3*(exp(ϵ_assoc[i,j][a,b]/T)-1)*κ[i,j][a,b]
# end

function  Δ(model::ePCSAFT, V, T, z,_data=@f(data))
    ϵ_assoc = model.params.epsilon_assoc.values
    κ = model.params.bondvol.values
    σ = model.params.sigma.values
    Δres = zero_assoc(κ,typeof(V+T+first(z)))
    for (idx,(i,j),(a,b)) in indices(Δres)
        if model.components[i]=="water" && model.components[j]=="water"
            σij = (2.792700+10.1100*exp(-0.01775*T)-1.41700*exp(-0.01146*T))*1e-10
        elseif model.components[i]=="water" && model.components[j]!="water"
            σij = ((2.792700+10.1100*exp(-0.01775*T)-1.41700*exp(-0.01146*T))*1e-10+model.params.sigma.values[j,j])/2
        elseif model.components[j]=="water" && model.components[i]!="water"
            σij = ((2.792700+10.1100*exp(-0.01775*T)-1.41700*exp(-0.01146*T))*1e-10+model.params.sigma.values[i,i])/2
        else
            σij = model.params.sigma.values[i,j]
        end
        gij = @f(g_hs,i,j,_data)
        Δres[idx] = gij*σij^3*(exp(ϵ_assoc[i,j][a,b]/T)-1)*κ[i,j][a,b]
    end
    return Δres
end