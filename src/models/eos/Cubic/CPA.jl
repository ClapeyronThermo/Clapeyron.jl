function a_res(model::CPAFamily,z,v,T)
    return a_SRK(model,z,v,T)+a_assoc(model,z,v,T)
end

function a_SRK(model::CPAFamily,z,v,T)
    x = z/sum(z[i] for i in model.components)
    n = sum(z)

    āᾱ = sum(sum(model.params.a[union(i,j)]*√(α(model,T,i)*α(model,T,j))*x[i]*x[j] for j in model.components) for i in model.components)
    b̄  = sum(sum(model.params.b[union(i,j)]*x[i]*x[j] for j in model.components) for i in model.components)
    return -log(1-n*b̄/v)-āᾱ/(R̄*T*b̄)*log(1+n*b̄/v)
end

function α(model::CPAFamily,T,i)
    Tc = model.params.Tc[i]
    c1  = model.params.c1[i]
    return (1+c1*(1-√(T/Tc)))^2
end

function a_assoc(model::CPAFamily, z, v, T)
    x = z/sum(z[i] for i in model.components)
    n_sites = model.params.n_sites

    X_iA = X_assoc(model,z,v,T)
    return sum(x[i]*sum(n_sites[i][a]*(log(X_iA[i,a])+(1-X_iA[i,a])/2) for a in keys(model.params.n_sites[i])) for i in model.components)
end

function X_assoc(model::CPAFamily, z, v, T)
    x = z/sum(z[i] for i in model.components)
    ρ = sum(z[i] for i in model.components)/v
    n_sites = model.params.n_sites
    X_iA = Dict()
    X_iA_old = Dict()
    tol = 1.
    iter = 1
    while tol > 1e-12 && iter < 100
        for i in model.components
            for a in keys(model.params.n_sites[i])
                A = 0.
                for j in model.components
                    B = 0
                    for b in keys(model.params.n_sites[j])
                        if haskey(model.params.epsilon_assoc,Set([(i,a),(j,b)]))
                            if iter!=1
                                B+=n_sites[j][b]*X_iA_old[j,b]*Δ(model,z,v,T,i,j,a,b)
                            else
                                B+=n_sites[j][b]*Δ(model,z,v,T,i,j,a,b)
                            end
                        end
                    end
                    A += ρ*x[j]*B
                end
                if iter == 1
                    X_iA[i,a] =0.5+0.5*(1+A)^-1
                else
                    X_iA[i,a] =0.5*X_iA_old[i,a]+0.5*(1+A)^-1
                end
            end
        end
        if iter == 1
            tol = sqrt(sum(sum((1. -X_iA[i,a])^2 for a in keys(model.params.n_sites[i])) for i in model.components))
        else
            tol = sqrt(sum(sum((X_iA_old[i,a] -X_iA[i,a])^2 for a in keys(model.params.n_sites[i])) for i in model.components))
        end
        X_iA_old = deepcopy(X_iA)
        iter += 1
    end

    return X_iA
end

function Δ(model::CPAFamily, z, v, T, i, j, a, b)
    ϵ_assoc = model.params.epsilon_assoc[Set([(i,a),(j,b)])]*1e2/R̄
    β = model.params.bond_vol[Set([(i,a),(j,b)])]*1e-3
    x = z/sum(z[i] for i in model.components)
    b̄  = sum(sum(model.params.b[union(i,j)]*x[i]*x[j] for j in model.components) for i in model.components)
    η = b̄*sum(z[i] for i in model.components)/(4*v)
    g = (1-1.9η)^-1
    b_ij = (model.params.b[i]+model.params.b[j])/2
    return g*(exp(ϵ_assoc/T)-1)*β*b_ij
end
