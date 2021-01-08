function a_res(model::sCKSAFTFamily, z, v, T)
    return a_seg(model,z,v,T) + a_chain(model,z,v,T) + a_assoc(model,z,v,T)
end

function a_seg(model::sCKSAFTFamily, z, v, T)
    x = z/sum(z[i] for i in model.components)
    m = model.params.segment
    m̄ = sum(x[i]*x[j]*(m[i]+m[j])/2 for i in model.components for j in model.components)
    return m̄*(a_hs(model,z,v,T) + a_disp(model,z,v,T))
end

function a_hs(model::sCKSAFTFamily, z, v, T)
    ζ0 = ζn(model,z,v,T, 0)
    ζ1 = ζn(model,z,v,T, 1)
    ζ2 = ζn(model,z,v,T, 2)
    ζ3 = ζn(model,z,v,T, 3)
    return 1/ζ0 * (3ζ1*ζ2/(1-ζ3) + ζ2^3/(ζ3*(1-ζ3)^2) + (ζ2^3/ζ3^2-ζ0)*log(1-ζ3))
end

function a_disp(model::sCKSAFTFamily, z, v, T)
    components = model.components
    x = z/sum(z[i] for i in components)
    m = model.params.segment
    m̄ = sum(x[i]*m[i] for i in components)
    vs = v/(N_A*sum(z)*m̄)
    v̄Ȳ = sum(x[i]*x[j]*m[i]*m[j]*(d(model,z,v,T,i,j)^3/√(2))*(exp(u(model,z,v,T,i,j)/T/2)-1) for i in components for j in components)/sum(x[i]*x[j]*m[i]*m[j] for i in components for j in components)
    return 36*log(vs/(vs+v̄Ȳ))
end

function d(model::sCKSAFTFamily, z, v, T, component)
    ϵ = model.params.epsilon[component]
    σ = model.params.sigma[component]
    return σ * (1 - 0.333exp(-3ϵ/T))
end

function d(model::sCKSAFTFamily, z, v, T, i,j)
    return (d(model::sCKSAFTFamily, z, v, T, i)+d(model::sCKSAFTFamily, z, v, T, j))/2
end

function u(model::sCKSAFTFamily, z, v, T, i,j)
    ϵ0 = model.params.epsilon[union(i,j)]
    return ϵ0*(1-10/T)
end

function ζn(model::sCKSAFTFamily, z, v, T, n)
    x = z/sum(z[i] for i in model.components)
    m = model.params.segment
    return N_A*sum(z[i] for i in model.components)*π/6/v * sum(x[i]*m[i]*d(model,z,v,T, i)^n for i in model.components)
end

function a_chain(model::sCKSAFTFamily, z, v, T)
    x = z/sum(z[i] for i in model.components)
    m = model.params.segment
    return sum(x[i]*(1-m[i])*log(g_hsij(model,z,v,T,i,i)) for i in model.components)
end

function g_hsij(model::sCKSAFTFamily, z, v, T, i, j)
    di = d(model,z,v,T, i)
    dj = d(model,z,v,T, j)
    ζ2 = ζn(model,z,v,T, 2)
    ζ3 = ζn(model,z,v,T, 3)
    return 1/(1-ζ3) + di*dj/(di+dj)*3ζ2/(1-ζ3)^2 + (di*dj/(di+dj))^2*2ζ2^2/(1-ζ3)^3
end

#constants should be allocated outside the function
const sCKSAFTconsts = (
    D =
    [0.9105631445 -0.3084016918 -0.0906148351;
    0.6361281449 0.1860531159 0.4527842806;
    2.6861347891 -2.5030047259 0.5962700728;
    -26.547362491 21.419793629 -1.7241829131;
    97.759208784 -65.255885330 -4.1302112531;
    -159.59154087 83.318680481 13.776631870;
    91.297774084 -33.746922930 -8.6728470368]
)


function a_assoc(model::sCKSAFTFamily, z, v, T)
    x = z/sum(z[i] for i in model.components)
    X_iA = X_assoc(model,z,v,T)
    return sum(x[i]*sum(log(X_iA[i,a])-X_iA[i,a]/2 + model.params.n_sites[i][a]/2 for a in keys(model.params.n_sites[i])) for i in model.components)
end

function X_assoc(model::sCKSAFTFamily, z, v, T)
    x = z/sum(z[i] for i in model.components)
    ρ = N_A*sum(z[i] for i in model.components)/v
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
                                B+=X_iA_old[j,b]*Δ(model,z,v,T,i,j,a,b)
                            else
                                B+=Δ(model,z,v,T,i,j,a,b)
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

function Δ(model::sCKSAFTFamily, z, v, T, i, j, a, b)
    ϵ_assoc = model.params.epsilon_assoc[Set([(i,a),(j,b)])]
    κ = model.params.bond_vol[Set([(i,a),(j,b)])]
    g = g_hsij(model,z,v,T,i,j)
    return g*d(model,z,v,T,i,j)^3*(exp(ϵ_assoc/T)-1)*κ
end
