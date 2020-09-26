const N_A = 6.02214086e23
const k_B = 1.38064852e-23
const R   = N_A*k_B

function a_res(model::PCSAFTFamily, z, v, T)
    return a_hc(model,z,v,T) + a_disp(model,z,v,T) #+ a_assoc(model,z,v,T)
end

function a_hc(model::PCSAFTFamily, z, v, T)
    x = z/sum(z[i] for i in model.components)
    m = model.parameters.segment
    m̄ = sum(x[i]*m[i] for i in model.components)
    return m̄*a_hs(model,z,v,T) - sum(x[i]*(m[i]-1)*log(g_hsij(model,z,v,T, i, i)) for i in model.components)
end

function a_disp(model::PCSAFTFamily, z, v, T)
    x = z/sum(z[i] for i in model.components)
    m = model.parameters.segment
    m̄ = sum(x[i]*m[i] for i in model.components)
    return -2*π*N_A*sum(z[i] for i in model.components)/v*I_n(model,z,v,T, 1)*m2ϵσ3(model,z,v,T, 1) - π*m̄*N_A*sum(z[i] for i in model.components)/v*C1(model,z,v,T)*I_n(model,z,v,T, 2)*m2ϵσ3(model,z,v,T, 2)
end

function d(model::PCSAFTFamily, z, v, T, component)
    ϵ = model.parameters.epsilon[component]
    σ = model.parameters.sigma[component]
    return σ * (1 - 0.12exp(-3ϵ/T))
end

function ζn(model::PCSAFTFamily, z, v, T, n)
    x = z/sum(z[i] for i in model.components)
    m = model.parameters.segment
    return N_A*sum(z[i] for i in model.components)*π/6/v * sum(x[i]*m[i]*d(model,z,v,T, i)^n for i in model.components)
end

function g_hsij(model::PCSAFTFamily, z, v, T, i, j)
    di = d(model,z,v,T, i)
    dj = d(model,z,v,T, j)
    ζ2 = ζn(model,z,v,T, 2)
    ζ3 = ζn(model,z,v,T, 3)
    return 1/(1-ζ3) + di*dj/(di+dj)*3ζ2/(1-ζ3)^2 + (di*dj/(di+dj))^2*2ζ2^2/(1-ζ3)^3
end

function a_hs(model::PCSAFTFamily, z, v, T)
    ζ0 = ζn(model,z,v,T, 0)
    ζ1 = ζn(model,z,v,T, 1)
    ζ2 = ζn(model,z,v,T, 2)
    ζ3 = ζn(model,z,v,T, 3)
    return 1/ζ0 * (3ζ1*ζ2/(1-ζ3) + ζ2^3/(ζ3*(1-ζ3)^2) + (ζ2^3/ζ3^2-ζ0)*log(1-ζ3))
end

function C1(model::PCSAFTFamily, z, v, T)
    x = z/sum(z[i] for i in model.components)
    η = ζn(model,z,v,T, 3)
    m = model.parameters.segment
    m̄ = sum(x[i]*m[i] for i in model.components)
    return (1 + m̄*(8η-2η^2)/(1-η)^4 + (1-m̄)*(20η-27η^2+12η^3-2η^4)/((1-η)*(2-η))^2)^-1
end

function m2ϵσ3(model::PCSAFTFamily, z, v, T, ϵ_power = 1)
    x = z/sum(z[i] for i in model.components)
    m = model.parameters.segment
    σ = model.parameters.sigma
    ϵ = model.parameters.epsilon
    k = model.parameters.k
    #= return sum(x[i]*x[j]*m[i]*m[j] * (sqrt(ϵ[i]*ϵ[j])*(1-k[union(i,j)])/T)^ϵ_power * (0.5*(σ[i]+σ[j]))^3 for i in model.components, j in model.components) =#
    return sum(x[i]*x[j]*m[i]*m[j] * (sqrt(ϵ[i]*ϵ[j])*(1)/T)^ϵ_power * (0.5*(σ[i]+σ[j]))^3 for i in model.components, j in model.components)
end

function I_n(model::PCSAFTFamily, z, v, T, n)
    x = z/sum(z[i] for i in model.components)
    m = model.parameters.segment
    m̄ = sum(x[i]*m[i] for i in model.components)
    η = ζn(model,z,v,T, 3)
    if n == 1
        corr = [0.9105631445 -0.3084016918 -0.0906148351;
                0.6361281449 0.1860531159 0.4527842806;
                2.6861347891 -2.5030047259 0.5962700728;
                -26.547362491 21.419793629 -1.7241829131;
                97.759208784 -65.255885330 -4.1302112531;
                -159.59154087 83.318680481 13.776631870;
                91.297774084 -33.746922930 -8.6728470368]
    elseif n == 2
        corr = [0.7240946941 -0.5755498075 0.0976883116;
                2.2382791861 0.6995095521 -0.2557574982;
                -4.0025849485 3.8925673390 -9.1558561530;
                -21.003576815 -17.215471648 20.642075974;
                26.855641363 192.67226447 -38.804430052;
                206.55133841 -161.82646165 93.626774077;
                -355.60235612 -165.20769346 -29.666905585]
    end
    return sum((corr[i+1,1] + (m̄-1)/m̄*corr[i+1,2] + (m̄-1)/m̄*(m̄-2)/m̄*corr[i+1,3]) * η^i for i = 0:6)
end

function a_assoc(model::PCSAFTFamily, z, v, T)
    x = z/sum(z[i] for i in model.components)
    M = model.sites
    X_iA = X_assoc(model,z,v,T)
    return sum(x[i]*sum(log(X_iA[i,a])-X_iA[i,a]/2 for a in 1:length(model.sites))+length(model.sites)/2 for i in model.components)
end

function X_assoc(model, z, v, T)
    x = z/sum(z[i] for i in model.components)
    ρ = N_A*sum(z[i] for i in model.components)/v
    X_iA = Matrix(undef,length(model.components),length(model.sites))
    X_iA_old = ones(length(model.components),length(model.sites))
    tol = 1.
    while tol > 1e-12
        for i in model.components
            for a in 1:length(model.sites)
                X_iA[i,a] =0.5*X_iA_old[i,a]+0.5*(1+sum(ρ*x[j]*sum(X_iA_old[j,b]*Δ(model,z,v,T,i,j,a,b) for b in 1:length(model.sites)) for j in model.components))^-1
            end
        end
        tol = norm(X_iA_old-X_iA)
        X_iA_old = deepcopy(X_iA)
    end
    return X_iA
end

function Δ(model, z, v, T, i, j, a, b)
    ϵ_assoc = model.parameters.epsilon_assoc[i,j,a,b]
    κ = model.parameters.Bond_vol[i,j,a,b]
    σ = model.parameters.sigma[i,j]
    g = g_hsij(model,z,v,T,i,j)
    return g*σ^3*(exp(ϵ_assoc/T)-1)*κ
end
