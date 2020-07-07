NParticles = 6.022140857e23;

function a_ideal(model::SAFTVRMie, conditions)
    0
end
function a_mono(model::SAFTVRMie, conditions)
    return a_hs(model,conditions)+a_disp(model,conditions)
end
function a_disp(model::SAFTVRMie, conditions)
    return a_1(model,conditions)+a_2(model,conditions)+a_3(model,conditions)
end

function a_hs(model::SAFTVRMie, conditions)
    ζ0   = ζn(model, conditions, 0)
    ζ1   = ζn(model, conditions, 1)
    ζ2   = ζn(model, conditions, 2)
    ζ3   = ζn(model, conditions, 3)
    Vol  = conditions.volume
    return 6/π/ρs(model,conditions)*(3ζ1*ζ2/(1-ζ3) + ζ2^3/(ζ3*(1-ζ3)^2) + (ζ2^3/ζ3^2-ζ0)*log(1-ζ3))
end

function ζn(model::SAFTVRMie, conditions,n)
    return π/6*ρs(model,conditions)*sum(xS(model,conditions,i)*d(model,conditions,i)^n for i in model.components)
end

function ρs(model::SAFTVRMie, conditions)
    Vol = conditions.volume
    x   = conditions.components
    m   = model.parameters.segment
    m̄   = sum(x[i]*m[i] for i in model.components)
    return NParticles/Vol*m̄
end

function xS(model::SAFTVRMie, conditions,i)
    x   = conditions.components
    m   = model.parameters.segment
    m̄   = sum(x[j]*m[j] for j in model.components)
    return x[i]*m[i]/m̄
end

function d(model::SAFTVRMie, conditions,i)
    Temp        = conditions.temperature
    ϵ           = model.parameters.epsilon[(i,i)]
    σ           = model.parameters.sigma[(i,i)]
    λR          = model.parameters.lambdaR[(i,i)]
    λA          = model.parameters.lambdaA[(i,i)]
    u           = [0.26356031971814109102031,1.41340305910651679221800,3.59642577104072208122300,7.08581000585883755692200,12.6408008442757826594300]
    w           = [0.5217556105828086524759,0.3986668110831759274500,7.5942449681707595390e-2,3.6117586799220484545e-3,2.3369972385776227891e-5]
    θ           = CMie(model,conditions,i,i)*ϵ/Temp
    return σ*(1-sum(w[j]*(θ./(θ+u[j]))^(1/λR)*(exp(θ*(1/(θ./(θ+u[j]))^(λA/λR)-1))/(u[j]+θ)/λR) for j in 1:5))
end

function CMie(model::SAFTVRMie, conditions,i,j)
    λR          = model.parameters.lambdaR
    λA          = model.parameters.lambdaA
    return (λR[(i,j)]/(λR[(i,j)]-λA[(i,j)]))*(λR[(i,j)]/λA[(i,j)])^(λA[(i,j)]/(λR[(i,j)]-λA[(i,j)]))
end

function ζ_x(model::SAFTVRMie, conditions)
    return pi/6*ρs(model,conditions)*sum(xS(model,conditions,i)*xS(model,conditions,j)*(d(model,conditions,i)+d(model,conditions,j))^3/8 for i in model.components for j in model.components)
end

function a_1(model::SAFTVRMie, conditions)
    x    = conditions.components
    Temp = conditions.temperature
    m    = model.parameters.segment
    m̄    = sum(x[j]*m[j] for j in model.components)
    return m̄/Temp*sum(xS(model,conditions,i)*xS(model,conditions,j)*a_1ij(model,conditions,i,j) for i in model.components for j in model.components)
end

function a_1ij(model::SAFTVRMie, conditions,i,j)
    ϵ           = model.parameters.epsilon[(i,j)]
    λR          = model.parameters.lambdaR[(i,j)]
    λA          = model.parameters.lambdaA[(i,j)]

    HSd         = (d(model,conditions,i)+d(model,conditions,j))/2
    x_0         = x0(model,conditions,i,j)
    C           = CMie(model,conditions,i,j)
    return 2*π*ϵ*HSd^3*C*ρs(model,conditions)*(x_0^λA*(a1s(model,conditions,λA)+B(model,conditions,λA,x_0))-x_0^λR*(a1s(model,conditions,λR)+B(model,conditions,λR,x_0)))
end

function a1s(model::SAFTVRMie, conditions,λ)
    ζeff = ζ_eff(model,conditions,λ)
    return -1/(λ-3)*(1-ζeff/2)/(1-ζeff)^3
end

function ζ_eff(model::SAFTVRMie, conditions,λ)
    ζx = ζ_x(model,conditions)
    return sum(c(model,conditions,λ,n)*ζx^n for n in 1:4)
end

function c(model::SAFTVRMie, conditions,λ,n)
    A           = [[0.81096   1.7888  -37.578   92.284],
                   [1.02050  -19.341   151.26  -463.50],
                  [-1.90570   22.845  -228.14   973.92],
                   [1.08850  -6.1962   106.98  -677.64]]
    return sum(A[n][m]/λ^(m-1) for m in 1:4)
end

function B(model::SAFTVRMie, conditions,λ,x0)
    I  = (1-x0^(3-λ))/(λ-3)
    J  = (1-(λ-3)*x0^(4-λ)+(λ-4)*x0^(3-λ))/((λ-3)*(λ-4))
    ζx = ζ_x(model,conditions)
    return I*(1-ζx/2)/(1-ζx)^3-9*J*ζx*(ζx+1)/(2*(1-ζx)^3)
end

function x0(model::SAFTVRMie, conditions,i,j)
    σ           = model.parameters.sigma
    HSd         = (d(model,conditions,i)+d(model,conditions,j))/2
    return σ[(i,j)]/HSd
end

function a_2(model::SAFTVRMie, conditions)
    x    = conditions.components
    Temp = conditions.temperature
    m    = model.parameters.segment
    m̄    = sum(x[j]*m[j] for j in model.components)
    return m̄/Temp^2*sum(xS(model,conditions,i)*xS(model,conditions,j)*a_2ij(model,conditions,i,j) for i in model.components for j in model.components)
end

function a_2ij(model::SAFTVRMie, conditions,i,j)
    ϵ           = model.parameters.epsilon[(i,j)]
    λR          = model.parameters.lambdaR[(i,j)]
    λA          = model.parameters.lambdaA[(i,j)]

    HSd         = (d(model,conditions,i)+d(model,conditions,j))/2
    x_0         = x0(model,conditions,i,j)
    C           = CMie(model,conditions,i,j)
    ζx          = ζ_x(model,conditions)
    KHS         = (1-ζx)^4/(1+4*ζx+4*ζx^2-4*ζx^3+ζx^4)
    return π*KHS*(1+χ(model,conditions,i,j))*ρs(model,conditions)*ϵ^2*HSd^3*C^2*(x_0^(2*λA)*(a1s(model,conditions,2*λA)+B(model,conditions,2*λA,x_0))-2*x_0^(λA+λR)*(a1s(model,conditions,λA+λR)+B(model,conditions,λA+λR,x_0))+x_0^(2*λR)*(a1s(model,conditions,2*λR)+B(model,conditions,2*λR,x_0)))
end

function χ(model::SAFTVRMie, conditions,i,j)
    λR          = model.parameters.lambdaR[(i,j)]
    λA          = model.parameters.lambdaA[(i,j)]

    ζstar = ζ_star(model, conditions)
    α     = CMie(model,conditions,i,j)*(1/(λA-3)-1/(λR-3))
    return f(α,1)*ζstar+f(α,2)*ζstar^5+f(α,3)*ζstar^8
end

function f(α,m)
    phi         = [[7.5365557, -359.440,  1550.9, -1.199320, -1911.2800,  9236.9,  10.0],
                  [-37.604630,  1825.60, -5070.1,  9.063632,  21390.175, -129430,  10.0],
                   [71.745953, -3168.00,  6534.6, -17.94820, -51320.700,  357230,  0.57],
                  [-46.835520,  1884.20, -3288.7,  11.34027,  37064.540, -315530, -6.70],
                  [-2.4679820,- 0.82376, -2.7171,  20.52142,  1103.7420,  1390.2, -8.00],
                  [-0.5027200, -3.19350,  2.0883, -56.63770, -3264.6100, -4518.2,   NaN],
                  [8.0956883,  3.70900,  0.0000,  40.53683,  2556.1810,  4241.6,   NaN]];
    return sum(phi[i+1][m]*α^i for i in 0:3)/(1+sum(phi[i+1][m]*α^(i-3) for i in 4:6))
end

function ζ_star(model::SAFTVRMie, conditions)
    σ  = model.parameters.sigma
    return ρs(model,conditions)*π/6*sum(xS(model,conditions,i)*xS(model,conditions,j)*(σ[(i,j)])^3 for i in model.components for j in model.components)
end

function a_3(model::SAFTVRMie, conditions)
    x    = conditions.components
    Temp = conditions.temperature
    m    = model.parameters.segment
    m̄    = sum(x[j]*m[j] for j in model.components)
    return m̄/Temp^3*sum(xS(model,conditions,i)*xS(model,conditions,j)*a_3ij(model,conditions,i,j) for i in model.components for j in model.components)
end

function a_3ij(model::SAFTVRMie, conditions,i,j)
    ϵ     = model.parameters.epsilon[(i,j)]
    λR    = model.parameters.lambdaR[(i,j)]
    λA    = model.parameters.lambdaA[(i,j)]
    ζstar = ζ_star(model, conditions)
    α     = CMie(model,conditions,i,j)*(1/(λA-3)-1/(λR-3))
    return -ϵ^3*f(α,4)*ζstar*exp(f(α,5)*ζstar+f(α,6)*ζstar^2)
end

function a_chain(model::SAFTVRMie, conditions)
    0
end
function a_assoc(model::SAFTVRMie, conditions)
    0
end
