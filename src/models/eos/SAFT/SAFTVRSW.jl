

function a_res(model::SAFTVRSWFamily, z,Vol,Temp)
    return a_mono(model,z,Vol,Temp)+a_chain(model,z,Vol,Temp)+a_assoc(model,z,Vol,Temp)
end

function a_mono(model::SAFTVRSWFamily, z,Vol,Temp)
    return a_hs(model,z,Vol,Temp)+a_disp(model,z,Vol,Temp)
end
function a_disp(model::SAFTVRSWFamily, z,Vol,Temp)
    return a_1(model,z,Vol,Temp)+a_2(model,z,Vol,Temp)
end

function a_hs(model::SAFTVRSWFamily, z,Vol,Temp)
    ζ0   = ζn(model, z,Vol,Temp, 0)
    ζ1   = ζn(model, z,Vol,Temp, 1)
    ζ2   = ζn(model, z,Vol,Temp, 2)
    ζ3   = ζn(model, z,Vol,Temp, 3)
    NParticles = N_A*sum(z[i] for i in model.components)
    x   = z/sum(z[i] for i in model.components)
    m   = model.params.segment
    m̄   = sum(x[i]*m[i] for i in model.components)
    return m̄*6/π/ρs(model,z,Vol,Temp)*(3ζ1*ζ2/(1-ζ3) + ζ2^3/(ζ3*(1-ζ3)^2) + (ζ2^3/ζ3^2-ζ0)*log(1-ζ3))
end

function ζn(model::SAFTVRSWFamily, z,Vol,Temp,n)
    σ = model.params.sigma
    return π/6*ρs(model,z,Vol,Temp)*sum(xS(model,z,Vol,Temp,i)*σ[i]^n for i in model.components)
end

function ρs(model::SAFTVRSWFamily, z,Vol,Temp)
    NParticles = N_A*sum(z[i] for i in model.components)
    x   = z/sum(z[i] for i in model.components)
    m   = model.params.segment
    m̄   = sum(x[i]*m[i] for i in model.components)
    return NParticles/Vol*m̄
end

function xS(model::SAFTVRSWFamily, z,Vol,Temp,i)
    x   = z/sum(z[i] for i in model.components)
    m   = model.params.segment
    m̄   = sum(x[j]*m[j] for j in model.components)
    return x[i]*m[i]/m̄
end

function ζ_x(model::SAFTVRSWFamily, z,Vol,Temp)
    σ = model.params.sigma
    return pi/6*ρs(model,z,Vol,Temp)*sum(xS(model,z,Vol,Temp,i)*xS(model,z,Vol,Temp,j)*(σ[union(i,j)])^3 for i in model.components for j in model.components)
end

function a_1(model::SAFTVRSWFamily, z,Vol,Temp)
    x    = z/sum(z[i] for i in model.components)

    m    = model.params.segment
    m̄    = sum(x[j]*m[j] for j in model.components)
    return -m̄/Temp*ρs(model,z,Vol,Temp)*sum(xS(model,z,Vol,Temp,i)*xS(model,z,Vol,Temp,j)*a_1ij(model,z,Vol,Temp,i,j) for i in model.components for j in model.components)
end

function a_1ij(model::SAFTVRSWFamily, z,Vol,Temp,i,j)
    ϵ           = model.params.epsilon[union(i,j)]
    λ           = model.params.lambda[union(i,j)]
    σ           = model.params.sigma[union(i,j)]

    α           = 2π*ϵ*σ^3*(λ^3-1)/3
    gHS0        = g_HS0(model,z,Vol,Temp,i,j)
    return α*gHS0
end

function ζ_eff(model::SAFTVRSWFamily, z,Vol,Temp,λ)
    ζx = ζ_x(model,z,Vol,Temp)
    return sum(c(model,z,Vol,Temp,λ,n)*ζx^n for n in 1:3)
end

function c(model::SAFTVRSWFamily, z,Vol,Temp,λ,n)
    A = [[2.25855   -1.50349  0.249434],
    [-0.66927  1.40049   -0.827739],
    [10.1576   -15.0427   5.30827]]
    return sum(A[n][m]*λ^(m-1) for m in 1:3)
end

function g_HS0(model::SAFTVRSWFamily,z,Vol,Temp,i,j)
    λ = model.params.lambda[union(i,j)]
    ζeff = ζ_eff(model,z,Vol,Temp,λ)
    return (1-ζeff/2)/(1-ζeff)^3
end

function a_2(model::SAFTVRSWFamily, z,Vol,Temp)
    x    = z/sum(z[i] for i in model.components)

    m    = model.params.segment
    m̄    = sum(x[j]*m[j] for j in model.components)
    return m̄/Temp^2*sum(xS(model,z,Vol,Temp,i)*xS(model,z,Vol,Temp,j)*a_2ij(model,z,Vol,Temp,i,j) for i in model.components for j in model.components)
end

function a_2ij(model::SAFTVRSWFamily, z,Vol,Temp,i,j)
    ϵ           = model.params.epsilon[union(i,j)]

    ζ0          = ζn(model, z,Vol,Temp, 0)
    ζ1          = ζn(model, z,Vol,Temp, 1)
    ζ2          = ζn(model, z,Vol,Temp, 2)
    ζ3          = ζn(model, z,Vol,Temp, 3)
    KHS         = ζ0*(1-ζ3)^4/(ζ0*(1-ζ3)^2+6*ζ1*ζ2*(1-ζ3)+9*ζ2^3)
    return 1/2*KHS*ϵ*ρs(model,z,Vol,Temp)*∂a1∂ρs(model,z,Vol,Temp,i,j)
end

function ∂a1∂ρs(model::SAFTVRSWFamily, z, Vol, Temp, i, j)
    ϵ           = model.params.epsilon[union(i,j)]
    λ           = model.params.lambda[union(i,j)]
    σ           = model.params.sigma[union(i,j)]

    α           = 2π*ϵ*σ^3*(λ^3-1)/3
    gHS0        = g_HS0(model,z,Vol,Temp,i,j)
    ρ_s         = ρs(model,z,Vol,Temp)
    ζeff        = ζ_eff(model,z,Vol,Temp,λ)
    ζx          = ζ_x(model,z,Vol,Temp)
    ∂ζeff∂ζx    = sum(n*c(model,z,Vol,Temp,λ,n)*ζx^n for n in 1:3)
    return -α*(gHS0+(5/2-ζeff)/(1-ζeff)^4*∂ζeff∂ζx)
end

function a_chain(model::SAFTVRSWFamily, z,Vol,Temp)
    x       = z/sum(z[i] for i in model.components)
    m       = model.params.segment
    return -sum(x[i]*(log(y_SW(model,z,Vol,Temp,i))*(m[i]-1)) for i in model.components)
end

function y_SW(model::SAFTVRSWFamily,z,Vol,Temp,i)
    ϵ    = model.params.epsilon[i]
    gSW  = g_SW(model,z,Vol,Temp,i,i)
    return gSW*exp(-ϵ/Temp)
end

function g_SW(model::SAFTVRSWFamily,z,Vol,Temp,i,j)
    ϵ    = model.params.epsilon[union(i,j)]
    gHS  = g_HS(model,z,Vol,Temp,i,j)
    g1   = g_1(model,z,Vol,Temp,i,j)
    return gHS+ϵ/Temp*g1
end

function g_HS(model::SAFTVRSWFamily,z,Vol,Temp,i,j)
    σ   = model.params.sigma
    ζ3  = ζn(model,z,Vol,Temp,3)
    D   = σ[i]*σ[j]/(σ[i]+σ[j])*sum(xS(model,z,Vol,Temp,k)*σ[k]^2 for k in model.components)/sum(xS(model,z,Vol,Temp,k)*σ[k]^3 for k in model.components)
    return 1/(1-ζ3)+3*D*ζ3/(1-ζ3)^2+2*(D*ζ3)^2/(1-ζ3)^3
end

function g_1(model::SAFTVRSWFamily,z,Vol,Temp,i,j)
    λ     = model.params.lambda[union(i,j)]
    gHS0  = g_HS0(model,z,Vol,Temp,i,j)
    ζx = ζ_x(model,z,Vol,Temp)
    ζeff = ζ_eff(model,z,Vol,Temp,λ)
    ∂ζeff∂ζx = sum(n*c(model,z,Vol,Temp,λ,n)*ζx^(n-1) for n in 1:3)
    ∂ζeff∂λ  = sum(∂c∂λ(model,z,Vol,Temp,λ,n)*ζx^n for n in 1:3)
    return gHS0+(λ^3-1)*(5/2-ζeff)/(1-ζeff)^4*(λ/3*∂ζeff∂λ-ζx*∂ζeff∂ζx)
end

function ∂c∂λ(model::SAFTVRSWFamily, z,Vol,Temp,λ,n)
    A = [[2.25855   -1.50349  0.249434],
    [-0.66927  1.40049   -0.827739],
    [10.1576   -15.0427   5.30827]]
    return sum((m-1)*A[n][m]*λ^(m-2) for m in 2:3)
end

function a_assoc(model::SAFTVRSWFamily, z, v, T)
    x = z/sum(z[i] for i in model.components)
    n_sites = model.params.n_sites
    X_iA = X_assoc(model,z,v,T)
    return sum(x[i]*sum(n_sites[i][a]*(log(X_iA[i,a])+(1-X_iA[i,a])/2) for a in keys(model.params.n_sites[i])) for i in model.components)
end

function X_assoc(model::SAFTVRSWFamily, z, v, T)
    x = z/sum(z[i] for i in model.components)
    ρ = N_A*sum(z[i] for i in model.components)/v
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

function Δ(model::SAFTVRSWFamily, z, v, T, i, j, a, b)
    ϵ_assoc = model.params.epsilon_assoc[Set([(i,a),(j,b)])]
    κ = model.params.bond_vol[Set([(i,a),(j,b)])]
    g = g_SW(model,z,v,T,i,j)
    return g*(exp(ϵ_assoc/T)-1)*κ
end
