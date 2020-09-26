N_A = 6.022140857e23;
k_B = 1.38064852e-23
R   = N_A*k_B
function a_res(model::SAFTFamily, z,Vol,Temp)
    return a_mono(model,z,Vol,Temp)+a_chain(model,z,Vol,Temp)
end

function a_mono(model::SAFTFamily, z,Vol,Temp)
    return a_hs(model,z,Vol,Temp)+a_disp(model,z,Vol,Temp)
end
function a_disp(model::SAFTFamily, z,Vol,Temp)
    return a_1(model,z,Vol,Temp)+a_2(model,z,Vol,Temp)+a_3(model,z,Vol,Temp)
end

function a_hs(model::SAFTFamily, z,Vol,Temp)
    ζ0   = ζn(model, z,Vol,Temp, 0)
    ζ1   = ζn(model, z,Vol,Temp, 1)
    ζ2   = ζn(model, z,Vol,Temp, 2)
    ζ3   = ζn(model, z,Vol,Temp, 3)
    NParticles = N_A*sum(z[i] for i in model.components)
    return 6*Vol/π/NParticles*(3ζ1*ζ2/(1-ζ3) + ζ2^3/(ζ3*(1-ζ3)^2) + (ζ2^3/ζ3^2-ζ0)*log(1-ζ3))
end

function ζn(model::SAFTFamily, z,Vol,Temp,n)
    return π/6*ρs(model,z,Vol,Temp)*sum(xS(model,z,Vol,Temp,i)*d(model,z,Vol,Temp,i)^n for i in model.components)
end

function ρs(model::SAFTFamily, z,Vol,Temp)
    NParticles = N_A*sum(z[i] for i in model.components)
    x   = z/sum(z[i] for i in model.components)
    m   = model.parameters.segment
    m̄   = sum(x[i]*m[i] for i in model.components)
    return NParticles/Vol*m̄
end

function xS(model::SAFTFamily, z,Vol,Temp,i)
    x   = z/sum(z[i] for i in model.components)
    m   = model.parameters.segment
    m̄   = sum(x[j]*m[j] for j in model.components)
    return x[i]*m[i]/m̄
end

function d(model::SAFTFamily, z,Vol,Temp,i)
    Temp        = Temp
    ϵ           = model.parameters.epsilon[(i,i)]
    σ           = model.parameters.sigma[(i,i)]
    λR          = model.parameters.lambdaR[(i,i)]
    λA          = model.parameters.lambdaA[(i,i)]
    u           = [0.26356031971814109102031,1.41340305910651679221800,3.59642577104072208122300,7.08581000585883755692200,12.6408008442757826594300]
    w           = [0.5217556105828086524759,0.3986668110831759274500,7.5942449681707595390e-2,3.6117586799220484545e-3,2.3369972385776227891e-5]
    θ           = CMie(model,z,Vol,Temp,i,i)*ϵ/Temp
    return σ*(1-sum(w[j]*(θ./(θ+u[j]))^(1/λR)*(exp(θ*(1/(θ./(θ+u[j]))^(λA/λR)-1))/(u[j]+θ)/λR) for j in 1:5))
end

function CMie(model::SAFTFamily, z,Vol,Temp,i,j)
    λR          = model.parameters.lambdaR
    λA          = model.parameters.lambdaA
    return (λR[(i,j)]/(λR[(i,j)]-λA[(i,j)]))*(λR[(i,j)]/λA[(i,j)])^(λA[(i,j)]/(λR[(i,j)]-λA[(i,j)]))
end

function ζ_x(model::SAFTFamily, z,Vol,Temp)
    return pi/6*ρs(model,z,Vol,Temp)*sum(xS(model,z,Vol,Temp,i)*xS(model,z,Vol,Temp,j)*(d(model,z,Vol,Temp,i)+d(model,z,Vol,Temp,j))^3/8 for i in model.components for j in model.components)
end

function a_1(model::SAFTFamily, z,Vol,Temp)
    x    = z/sum(z[i] for i in model.components)

    m    = model.parameters.segment
    m̄    = sum(x[j]*m[j] for j in model.components)
    return m̄/Temp*sum(xS(model,z,Vol,Temp,i)*xS(model,z,Vol,Temp,j)*a_1ij(model,z,Vol,Temp,i,j) for i in model.components for j in model.components)
end

function a_1ij(model::SAFTFamily, z,Vol,Temp,i,j)
    ϵ           = model.parameters.epsilon[(i,j)]
    λR          = model.parameters.lambdaR[(i,j)]
    λA          = model.parameters.lambdaA[(i,j)]

    HSd         = (d(model,z,Vol,Temp,i)+d(model,z,Vol,Temp,j))/2
    x_0         = x0(model,z,Vol,Temp,i,j)
    C           = CMie(model,z,Vol,Temp,i,j)
    return 2*π*ϵ*HSd^3*C*ρs(model,z,Vol,Temp)*(x_0^λA*(a1s(model,z,Vol,Temp,λA)+B(model,z,Vol,Temp,λA,x_0))-x_0^λR*(a1s(model,z,Vol,Temp,λR)+B(model,z,Vol,Temp,λR,x_0)))
end

function a1s(model::SAFTFamily, z,Vol,Temp,λ)
    ζeff = ζ_eff(model,z,Vol,Temp,λ)
    return -1/(λ-3)*(1-ζeff/2)/(1-ζeff)^3
end

function ζ_eff(model::SAFTFamily, z,Vol,Temp,λ)
    ζx = ζ_x(model,z,Vol,Temp)
    return sum(c(model,z,Vol,Temp,λ,n)*ζx^n for n in 1:4)
end

function c(model::SAFTFamily, z,Vol,Temp,λ,n)
    A           = [[0.81096   1.7888  -37.578   92.284],
                   [1.02050  -19.341   151.26  -463.50],
                  [-1.90570   22.845  -228.14   973.92],
                   [1.08850  -6.1962   106.98  -677.64]]
    return sum(A[n][m]/λ^(m-1) for m in 1:4)
end

function B(model::SAFTFamily, z,Vol,Temp,λ,x0)
    I  = (1-x0^(3-λ))/(λ-3)
    J  = (1-(λ-3)*x0^(4-λ)+(λ-4)*x0^(3-λ))/((λ-3)*(λ-4))
    ζx = ζ_x(model,z,Vol,Temp)
    return I*(1-ζx/2)/(1-ζx)^3-9*J*ζx*(ζx+1)/(2*(1-ζx)^3)
end

function x0(model::SAFTFamily, z,Vol,Temp,i,j)
    σ           = model.parameters.sigma
    HSd         = (d(model,z,Vol,Temp,i)+d(model,z,Vol,Temp,j))/2
    return σ[(i,j)]/HSd
end

function a_2(model::SAFTFamily, z,Vol,Temp)
    x    = z/sum(z[i] for i in model.components)

    m    = model.parameters.segment
    m̄    = sum(x[j]*m[j] for j in model.components)
    return m̄/Temp^2*sum(xS(model,z,Vol,Temp,i)*xS(model,z,Vol,Temp,j)*a_2ij(model,z,Vol,Temp,i,j) for i in model.components for j in model.components)
end

function a_2ij(model::SAFTFamily, z,Vol,Temp,i,j)
    ϵ           = model.parameters.epsilon[(i,j)]
    λR          = model.parameters.lambdaR[(i,j)]
    λA          = model.parameters.lambdaA[(i,j)]

    HSd         = (d(model,z,Vol,Temp,i)+d(model,z,Vol,Temp,j))/2
    x_0         = x0(model,z,Vol,Temp,i,j)
    C           = CMie(model,z,Vol,Temp,i,j)
    ζx          = ζ_x(model,z,Vol,Temp)
    KHS         = (1-ζx)^4/(1+4*ζx+4*ζx^2-4*ζx^3+ζx^4)
    return π*KHS*(1+χ(model,z,Vol,Temp,i,j))*ρs(model,z,Vol,Temp)*ϵ^2*HSd^3*C^2*(x_0^(2*λA)*(a1s(model,z,Vol,Temp,2*λA)+B(model,z,Vol,Temp,2*λA,x_0))-2*x_0^(λA+λR)*(a1s(model,z,Vol,Temp,λA+λR)+B(model,z,Vol,Temp,λA+λR,x_0))+x_0^(2*λR)*(a1s(model,z,Vol,Temp,2*λR)+B(model,z,Vol,Temp,2*λR,x_0)))
end

function χ(model::SAFTFamily, z,Vol,Temp,i,j)
    λR          = model.parameters.lambdaR[(i,j)]
    λA          = model.parameters.lambdaA[(i,j)]

    ζstar = ζ_star(model, z,Vol,Temp)
    α     = CMie(model,z,Vol,Temp,i,j)*(1/(λA-3)-1/(λR-3))
    return f(α,1)*ζstar+f(α,2)*ζstar^5+f(α,3)*ζstar^8
end

function f(α,m)
    ϕ           = [[7.5365557, -359.440,  1550.9, -1.199320, -1911.2800,  9236.9,  10.0],
                  [-37.604630,  1825.60, -5070.1,  9.063632,  21390.175, -129430,  10.0],
                   [71.745953, -3168.00,  6534.6, -17.94820, -51320.700,  357230,  0.57],
                  [-46.835520,  1884.20, -3288.7,  11.34027,  37064.540, -315530, -6.70],
                  [-2.4679820,- 0.82376, -2.7171,  20.52142,  1103.7420,  1390.2, -8.00],
                  [-0.5027200, -3.19350,  2.0883, -56.63770, -3264.6100, -4518.2,   NaN],
                  [8.0956883,  3.70900,  0.0000,  40.53683,  2556.1810,  4241.6,   NaN]];
    return sum(ϕ[i+1][m]*α^i for i in 0:3)/(1+sum(ϕ[i+1][m]*α^(i-3) for i in 4:6))
end

function ζ_star(model::SAFTFamily, z,Vol,Temp)
    σ  = model.parameters.sigma
    return ρs(model,z,Vol,Temp)*π/6*sum(xS(model,z,Vol,Temp,i)*xS(model,z,Vol,Temp,j)*(σ[(i,j)])^3 for i in model.components for j in model.components)
end

function a_3(model::SAFTFamily, z,Vol,Temp)
    x    = z/sum(z[i] for i in model.components)

    m    = model.parameters.segment
    m̄    = sum(x[j]*m[j] for j in model.components)
    return m̄/Temp^3*sum(xS(model,z,Vol,Temp,i)*xS(model,z,Vol,Temp,j)*a_3ij(model,z,Vol,Temp,i,j) for i in model.components for j in model.components)
end

function a_3ij(model::SAFTFamily, z,Vol,Temp,i,j)
    ϵ     = model.parameters.epsilon[(i,j)]
    λR    = model.parameters.lambdaR[(i,j)]
    λA    = model.parameters.lambdaA[(i,j)]
    ζstar = ζ_star(model, z,Vol,Temp)
    α     = CMie(model,z,Vol,Temp,i,j)*(1/(λA-3)-1/(λR-3))
    return -ϵ^3*f(α,4)*ζstar*exp(f(α,5)*ζstar+f(α,6)*ζstar^2)
end

function a_chain(model::SAFTFamily, z,Vol,Temp)
    x       = z/sum(z[i] for i in model.components)
    m       = model.parameters.segment
    return -sum(x[i]*(log(g_Mie(model,z,Vol,Temp,i))*(m[i]-1)) for i in model.components)
end

function g_Mie(model::SAFTFamily,z,Vol,Temp,i)
    ϵ    = model.parameters.epsilon[(i,i)]
    gHS  = g_HS(model,z,Vol,Temp,i)
    g1   = g_1(model,z,Vol,Temp,i)
    g2   = g_2(model,z,Vol,Temp,i)
    gMie = gHS*exp(ϵ/Temp*g1/gHS+(ϵ/Temp)^2*g2/gHS);
    return gMie
end

function g_HS(model::SAFTFamily,z,Vol,Temp,i)
    ζx  = ζ_x(model,z,Vol,Temp)
    x_0 = x0(model,z,Vol,Temp,i,i)
    k0  = -log(1-ζx)+(42ζx-39ζx^2+9ζx^3-2ζx^4)/(6*(1-ζx)^3)
    k1  = (ζx^4+6*ζx^2-12*ζx)/(2*(1-ζx)^3)
    k2  = -3*ζx^2/(8*(1-ζx)^2)
    k3  = (-ζx^4+3*ζx^2+3*ζx)/(6*(1-ζx)^3)
    gHS = exp(k0+x_0*k1+x_0^2*k2+x_0^3*k3)
    return gHS
end

function g_1(model::SAFTFamily,z,Vol,Temp,i)
    λR          = model.parameters.lambdaR[(i,i)]
    λA          = model.parameters.lambdaA[(i,i)]
    x_0 = x0(model,z,Vol,Temp,i,i)
    g1 = 3*da1(model,z,Vol,Temp,i)-CMie(model,z,Vol,Temp,i,i)*(λA*x_0^λA*(a1s(model,z,Vol,Temp,λA)+B(model,z,Vol,Temp,λA,x_0))-λR*x_0^λR*(a1s(model,z,Vol,Temp,λR)+B(model,z,Vol,Temp,λR,x_0)))
    return g1
end

function da1(model::SAFTFamily,z,Vol,Temp,i)
    λR  = model.parameters.lambdaR[(i,i)]
    λA  = model.parameters.lambdaA[(i,i)]
    x_0 = x0(model,z,Vol,Temp,i,i)
    return CMie(model,z,Vol,Temp,i,i)*(x_0^λA*(da1s(model,z,Vol,Temp,λA)+dB(model,z,Vol,Temp,λA,x_0))-x_0^λR*(da1s(model,z,Vol,Temp,λR)+dB(model,z,Vol,Temp,λR,x_0)))
end

function da1s(model::SAFTFamily,z,Vol,Temp,λ)
    ζeff  = ζ_eff(model,z,Vol,Temp,λ)
    dζeff = dζ_eff(model,z,Vol,Temp,λ)
    return -1/(λ-3)*((1-ζeff/2)/(1-ζeff)^3+ρs(model,z,Vol,Temp)*((3*(1-ζeff/2)*(1-ζeff)^2-0.5*(1-ζeff)^3)/(1-ζeff)^6)*dζeff);
end

function dζ_eff(model::SAFTFamily,z,Vol,Temp,λ)
    ζx = ζ_x(model,z,Vol,Temp)
    return sum(n*c(model,z,Vol,Temp,λ,n)*ζx^(n-1) for n in 1:4)*ζx/ρs(model,z,Vol,Temp)
end

function dB(model::SAFTFamily,z,Vol,Temp,λ,x_0)
    I  = (1-x_0^(3-λ))/(λ-3)
    J  = (1-(λ-3)*x_0^(4-λ)+(λ-4)*x_0^(3-λ))/((λ-3)*(λ-4))
    ζx = ζ_x(model,z,Vol,Temp)
    return (((1-ζx/2)*I/(1-ζx)^3-9*ζx*(1+ζx)*J/(2*(1-ζx)^3))+ρs(model,z,Vol,Temp)*((3*(1-ζx/2)*(1-ζx)^2
            -0.5*(1-ζx)^3)*I/(1-ζx)^6-9*J*((1+2*ζx)*(1-ζx)^3+ζx*(1+ζx)*3*(1-ζx)^2)/(2*(1-ζx)^6))*ζx/ρs(model,z,Vol,Temp));
end

function g_2(model::SAFTFamily,z,Vol,Temp,i)
        γ    = γ_corr(model::SAFTFamily,z,Vol,Temp,i)
        gMCA = g_MCA(model::SAFTFamily,z,Vol,Temp,i)
    return (1+γ)*gMCA
end

function γ_corr(model::SAFTFamily,z,Vol,Temp,i)
    ϕ    = [[7.5365557, -359.440,  1550.9, -1.199320, -1911.2800,  9236.9,  10.0],
           [-37.604630,  1825.60, -5070.1,  9.063632,  21390.175, -129430,  10.0],
            [71.745953, -3168.00,  6534.6, -17.94820, -51320.700,  357230,  0.57],
           [-46.835520,  1884.20, -3288.7,  11.34027,  37064.540, -315530, -6.70],
           [-2.4679820,- 0.82376, -2.7171,  20.52142,  1103.7420,  1390.2, -8.00],
           [-0.5027200, -3.19350,  2.0883, -56.63770, -3264.6100, -4518.2,   NaN],
            [8.0956883,  3.70900,  0.0000,  40.53683,  2556.1810,  4241.6,   NaN]];
   ϵ     = model.parameters.epsilon[(i,i)]
   λR    = model.parameters.lambdaR[(i,i)]
   λA    = model.parameters.lambdaA[(i,i)]
   ζstar = ζ_star(model, z,Vol,Temp)
   α     = CMie(model,z,Vol,Temp,i,i)*(1/(λA-3)-1/(λR-3))
   return ϕ[1][7]*(1-tanh(ϕ[2][7]*(ϕ[3][7]-α)))*ζstar*(exp(ϵ/Temp)-1)*exp(ϕ[4][7]*ζstar+ϕ[5][7]*ζstar^2)
end

function g_MCA(model::SAFTFamily,z,Vol,Temp,i)
    λR  = model.parameters.lambdaR[(i,i)]
    λA  = model.parameters.lambdaA[(i,i)]
    x_0 = x0(model,z,Vol,Temp,i,i)
    ζx  = ζ_x(model,z,Vol,Temp)
    KHS = (1-ζx)^4/(1+4*ζx+4*ζx^2-4*ζx^3+ζx^4)
    return 3*da2(model,z,Vol,Temp,i)-KHS*CMie(model,z,Vol,Temp,i,i)^2*(
          λR*x_0^(2*λR)*(a1s(model,z,Vol,Temp,2*λR)+B(model,z,Vol,Temp,2*λR,x_0))-
          (λA+λR)*x_0^(λA+λR)*(a1s(model,z,Vol,Temp,λA+λR)+B(model,z,Vol,Temp,λA+λR,x_0))+
          λA*x_0^(2*λA)*(a1s(model,z,Vol,Temp,2*λA)+B(model,z,Vol,Temp,2*λA,x_0)))
end
function da2(model::SAFTFamily,z,Vol,Temp,i)
    λR   = model.parameters.lambdaR[(i,i)]
    λA   = model.parameters.lambdaA[(i,i)]
    x_0  = x0(model,z,Vol,Temp,i,i)
    ζx   = ζ_x(model,z,Vol,Temp)
    KHS  = (1-ζx)^4/(1+4*ζx+4*ζx^2-4*ζx^3+ζx^4)
    dKHS = -((4*(1-ζx)^3*(1+4*ζx+4*ζx^2-4*ζx^3+ζx^4)+(1-ζx)^4*(4+8*ζx-12*ζx^2+4*ζx^3))/(1+4*ζx+4*ζx^2-4*ζx^3+ζx^4)^2)*ζx/ρs(model,z,Vol,Temp);
    return 0.5*CMie(model,z,Vol,Temp,i,i)^2*(
           ρs(model,z,Vol,Temp)*dKHS*(x_0^(2*λA)*(a1s(model,z,Vol,Temp,2*λA)+B(model,z,Vol,Temp,2*λA,x_0))-
           2*x_0^(λA+λR)*(a1s(model,z,Vol,Temp,λA+λR)+B(model,z,Vol,Temp,λA+λR,x_0))+
           x_0^(2*λR)*(a1s(model,z,Vol,Temp,2*λR)+B(model,z,Vol,Temp,2*λR,x_0)))+
           KHS*(x_0^(2*λA)*(da1s(model,z,Vol,Temp,2*λA)+dB(model,z,Vol,Temp,2*λA,x_0))-
           2*x_0^(λA+λR)*(da1s(model,z,Vol,Temp,λA+λR)+dB(model,z,Vol,Temp,λA+λR,x_0))+
           x_0^(2*λR)*(da1s(model,z,Vol,Temp,2*λR)+dB(model,z,Vol,Temp,2*λR,x_0))))
end

function a_assoc(model::SAFTFamily, z,Vol,Temp)
    0
end
