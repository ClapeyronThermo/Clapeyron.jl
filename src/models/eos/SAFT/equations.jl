function a_assoc(model::Union{SAFTModel,CPAModel}, V, T, z)
    if iszero(length(model.sites.flattenedsites))
        return zero(V+T+first(z))
    end
    X_ = @f(X)
    n = model.sites.n_sites
    res =  ∑(z[i]*∑(n[i][a] * (log(X_[i][a]) - X_[i][a]/2 + 0.5) for a ∈ @sites(i)) for i ∈ @comps)/sum(z)
    return res
end