function a_res(model::SAFTVRQMieFamily, z,Vol,Temp)
    return a_mono(model,z,Vol,Temp)
end

function a_mono(model::SAFTVRQMieFamily, z,Vol,Temp)
    return a_hs(model,z,Vol,Temp)+a_disp(model,z,Vol,Temp)
end
function a_disp(model::SAFTVRQMieFamily, z,Vol,Temp)
    return a_1(model,z,Vol,Temp)+a_2(model,z,Vol,Temp)+a_3(model,z,Vol,Temp)
end

function a_hs(model::SAFTVRQMieFamily, z,Vol,Temp)
    ζ0   = ζn(model, z,Vol,Temp, 0)
    ζ1   = ζn(model, z,Vol,Temp, 1)
    ζ2   = ζn(model, z,Vol,Temp, 2)
    ζ3   = ζn(model, z,Vol,Temp, 3)
    NParticles = N_A*sum(z[i] for i in model.components)
    return 6*Vol/π/NParticles*(3ζ1*ζ2/(1-ζ3) + ζ2^3/(ζ3*(1-ζ3)^2) + (ζ2^3/ζ3^2-ζ0)*log(1-ζ3))
end

function ζn(model::SAFTVRQMieFamily, z,Vol,Temp,n)
    return π/6*ρs(model,z,Vol,Temp)*sum(xS(model,z,Vol,Temp,i)*d(model,z,Vol,Temp,i)^n for i in model.components)
end

function ρs(model::SAFTVRQMieFamily, z,Vol,Temp)
    NParticles = N_A*sum(z[i] for i in model.components)
    x   = z/sum(z[i] for i in model.components)
    m   = model.params.segment
    m̄   = sum(x[i]*m[i] for i in model.components)
    return NParticles/Vol*m̄
end

function xS(model::SAFTVRQMieFamily, z,Vol,Temp,i)
    x   = z/sum(z[i] for i in model.components)
    m   = model.params.segment
    m̄   = sum(x[j]*m[j] for j in model.components)
    return x[i]*m[i]/m̄
end

function d(model::SAFTVRQMieFamily, z,Vol,Temp,i)
    ϵ           = model.params.epsilon[union(i,i)]
    σ           = model.params.sigma[union(i,i)]
    λR          = model.params.lambdaR[union(i,i)]
    λA          = model.params.lambdaA[union(i,i)]
    Mr          = model.params.MolarMass[i]/N_A

    D           = ħ^2/(12*k_B*Temp*Mr*σ^2)
    β           = ϵ/Temp
    u(x)        = CMie(model, z,Vol,Temp,i,i)*((x^-λR-x^-λA)+
                                        D*(Q1(λR)*x^-(λR+2)-Q1(λA)*x^-(λA+2))+
                                        D^2*(Q2(λR)*x^-(λR+4)-Q2(λA)*x^-(λA+4)))
    du(x)       = -CMie(model, z,Vol,Temp,i,i)*((λR*x^-(λR+1)-λA*x^-(λA+1))+
                                        D*(Q1(λR)*(λR+2)*x^-(λR+3)-Q1(λA)*(λA+2)*x^-(λA+3))+
                                        D^2*(Q2(λR)*(λR+4)*x^-(λR+5)-Q2(λA)*(λA+4)*x^-(λA+5)))
    d2u(x)      = CMie(model, z,Vol,Temp,i,i)*((λR*(λR+1)*x^-(λR+2)-λA*(λA+1)*x^-(λA+2))+
                                        D*(Q1(λR)*(λR+2)*(λR+3)*x^-(λR+4)-Q1(λA)*(λA+2)*(λA+3)*x^-(λA+4))+
                                        D^2*(Q2(λR)*(λR+4)*(λR+5)*x^-(λR+6)-Q2(λA)*(λA+4)*(λA+5)*x^-(λA+6)))


    f(x)        = exp(-β*u(x))
    g(x)        = -β*exp(-β*u(x))*du(x)
    h(x)        = β*exp(-β*u(x))*(β*du(x)^2-d2u(x))
    x_min       = Halley(f,g,h,1.)
    σ_eff        = σeff(model::SAFTVRQMieFamily,z,Vol,Temp,i,i)/σ
    x           = [-0.9956571630258080807355,-0.973906528517171720078,-0.9301574913557082260012,-0.8650633666889845107321,-0.7808177265864168970637,-0.6794095682990244062343,-0.562757134668604683339,
                   -0.4333953941292471907993,-0.294392862701460198131,-0.1488743389816312108848,0,0.1488743389816312108848,0.2943928627014601981311,0.4333953941292471907993,
                    0.562757134668604683339,0.6794095682990244062343,0.7808177265864168970637,0.865063366688984510732,0.9301574913557082260012,0.973906528517171720078,0.9956571630258080807355]
    w           =[0.0116946388673718742781,0.0325581623079647274788,0.0547558965743519960314,0.075039674810919952767,0.093125454583697605535,0.1093871588022976418992,0.123491976262065851078,
                  0.134709217311473325928,0.142775938577060080797,0.1477391049013384913748,0.149445554002916905665,0.1477391049013384913748,0.1427759385770600807971,0.134709217311473325928,
                  0.123491976262065851078,0.109387158802297641899,0.093125454583697605535,0.075039674810919952767,0.05475589657435199603138,0.032558162307964727479,0.0116946388673718742781]
    return σ*(σ_eff-(σ_eff-x_min)/2*sum(w[i]*f((σ_eff-x_min)/2*x[i]+(σ_eff+x_min)/2) for i in 1:length(x)))
end

function σeff(model::SAFTVRQMieFamily,z,Vol,Temp,i,j)
    σ           = model.params.sigma[union(i,j)]
    λR          = model.params.lambdaR[union(i,j)]
    λA          = model.params.lambdaA[union(i,j)]
    Mr          = model.params.MolarMass[union(i,j)]/N_A

    D           = ħ^2/(12*k_B*Temp*Mr*σ^2)

    f(x) = D^2*Q2(λR)+D*Q1(λR)*x^2+x^4-x^(λR-λA)*(D^2*Q2(λA)+D*Q1(λA)*x^2+x^4)
    g(x) = 2*D*Q1(λR)*x+4x^3-x^(λR-λA-1)*((λR-λA)*D^2*Q2(λA)+(λR-λA+2)*D*Q1(λA)*x^2+(λR-λA+4)*x^4)
    h(x) = 2*D*Q1(λR)+12x^2-x^(λR-λA-2)*((λR-λA-1)*(λR-λA)*D^2*Q2(λA)+((λR-λA-1)*(λR-λA+2)+2*(λR-λA+2))*D*Q1(λA)*x^2+((λR-λA-1)*(λR-λA+4)+4*(λR-λA+4))*x^4)

    return σ*Halley(f,g,h,1.)
end

function ϵeff(model::SAFTVRQMieFamily,z,Vol,Temp,i,j)
    ϵ           = model.params.epsilon[union(i,j)]
    σ           = model.params.sigma[union(i,j)]
    λR          = model.params.lambdaR[union(i,j)]
    λA          = model.params.lambdaA[union(i,j)]
    Mr          = model.params.MolarMass[union(i,j)]/N_A

    D           = ħ^2/(12*k_B*Temp*Mr*σ^2)

    u(x) = CMie(model, z,Vol,Temp,i,j)*((x^-λR-x^-λA)+
                                        D*(Q1(λR)*x^-(λR+2)-Q1(λA)*x^-(λA+2))+
                                        D^2*(Q2(λR)*x^-(λR+4)-Q2(λA)*x^-(λA+4)))
    f(x) = (λR+4)*D^2*Q2(λR)+(λR+2)*D*Q1(λR)*x^2+λR*x^4-x^(λR-λA)*((λA+4)*D^2*Q2(λA)+(λA+2)*D*Q1(λA)*x^2+λA*x^4)
    g(x) = 2*(λR+2)*D*Q1(λR)*x+4λR*x^3-x^(λR-λA-1)*((λA+4)*(λR-λA)*D^2*Q2(λA)+(λA+2)*(λR-λA+2)*D*Q1(λA)*x^2+λA*(λR-λA+4)*x^4)
    h(x) = 2*(λR+2)*D*Q1(λR)+12λR*x^2-x^(λR-λA-2)*((λA+4)*(λR-λA-1)*(λR-λA)*D^2*Q2(λA)+(λA+2)*((λR-λA-1)*(λR-λA+2)+2*(λR-λA+2))*D*Q1(λA)*x^2+λA*((λR-λA-1)*(λR-λA+4)+4*(λR-λA+4))*x^4)

    σ_min = Halley(f,g,h,(λR/λA)^(1/(λR-λA)))
    return -ϵ*u(σ_min)
end

function Halley(f,g,h,x0)
    tolx = 1.
    tolf = 1.
    f0 = f(x0)
    if f0<1e-16
        return x0
    else
        while tolf>1e-16 && tolx>1e-15
            d = f(x0)/g(x0)*(1-f(x0)*h(x0)/(2*g(x0)^2))^-1
            x0=x0-d
            f0 = f(x0)
            tolf = abs(f0)
            tolx = abs(d)
        end
        return x0
    end
end

function CMie(model::SAFTVRQMieFamily, z,Vol,Temp,i,j)
    λR          = model.params.lambdaR
    λA          = model.params.lambdaA
    return (λR[union(i,j)]/(λR[union(i,j)]-λA[union(i,j)]))*(λR[union(i,j)]/λA[union(i,j)])^(λA[union(i,j)]/(λR[union(i,j)]-λA[union(i,j)]))
end

function ζ_x(model::SAFTVRQMieFamily, z,Vol,Temp)
    return pi/6*ρs(model,z,Vol,Temp)*sum(xS(model,z,Vol,Temp,i)*xS(model,z,Vol,Temp,j)*(d(model,z,Vol,Temp,i)+d(model,z,Vol,Temp,j))^3/8 for i in model.components for j in model.components)
end

function a_1(model::SAFTVRQMieFamily, z,Vol,Temp)
    x    = z/sum(z[i] for i in model.components)

    m    = model.params.segment
    m̄    = sum(x[j]*m[j] for j in model.components)
    return m̄/Temp*sum(xS(model,z,Vol,Temp,i)*xS(model,z,Vol,Temp,j)*a_1ij(model,z,Vol,Temp,i,j) for i in model.components for j in model.components)
end

function a_1ij(model::SAFTVRQMieFamily, z,Vol,Temp,i,j)
    ϵ           = model.params.epsilon[union(i,j)]
    σ           = model.params.sigma[union(i,j)]
    λR          = model.params.lambdaR[union(i,j)]
    λA          = model.params.lambdaA[union(i,j)]
    Mr          = model.params.MolarMass[union(i,j)]/N_A

    D           = ħ^2/(12*k_B*Temp*Mr*σ^2)
    HSd         = (d(model,z,Vol,Temp,i)+d(model,z,Vol,Temp,j))/2
    x_0eff      = x0eff(model,z,Vol,Temp,i,j)
    x_0         = x0(model,z,Vol,Temp,i,j)
    C           = CMie(model,z,Vol,Temp,i,j)
    return 2*π*ϵ*HSd^3*C*ρs(model,z,Vol,Temp)*((x_0^λA*(a1s(model,z,Vol,Temp,λA)+B(model,z,Vol,Temp,λA,x_0eff))-x_0^λR*(a1s(model,z,Vol,Temp,λR)+B(model,z,Vol,Temp,λR,x_0eff)))+
                                               (x_0^(λA+2)*Q1(λA)*(a1s(model,z,Vol,Temp,λA+2)+B(model,z,Vol,Temp,λA+2,x_0eff))-x_0^(λR+2)*Q1(λR)*(a1s(model,z,Vol,Temp,λR+2)+B(model,z,Vol,Temp,λR+2,x_0eff)))*D+
                                               (x_0^(λA+4)*Q2(λA)*(a1s(model,z,Vol,Temp,λA+4)+B(model,z,Vol,Temp,λA+4,x_0eff))-x_0^(λR+4)*Q2(λR)*(a1s(model,z,Vol,Temp,λR+4)+B(model,z,Vol,Temp,λR+4,x_0eff)))*D^2)
end

function Q1(λ)
    return λ*(λ-1)
end

function Q2(λ)
    return 1/2*(λ+2)*(λ+1)*λ*(λ-1)
end

function a1s(model::SAFTVRQMieFamily, z,Vol,Temp,λ)
    ζeff = ζ_eff(model,z,Vol,Temp,λ)
    return -1/(λ-3)*(1-ζeff/2)/(1-ζeff)^3
end

function ζ_eff(model::SAFTVRQMieFamily, z,Vol,Temp,λ)
    ζx = ζ_x(model,z,Vol,Temp)
    return sum(c(model,z,Vol,Temp,λ,n)*ζx^n for n in 1:4)
end

function c(model::SAFTVRQMieFamily, z,Vol,Temp,λ,n)
    A           = [[0.81096   1.7888  -37.578   92.284],
                   [1.02050  -19.341   151.26  -463.50],
                  [-1.90570   22.845  -228.14   973.92],
                   [1.08850  -6.1962   106.98  -677.64]]
    return sum(A[n][m]/λ^(m-1) for m in 1:4)
end

function B(model::SAFTVRQMieFamily, z,Vol,Temp,λ,x0)
    I  = (1-x0^(3-λ))/(λ-3)
    J  = (1-(λ-3)*x0^(4-λ)+(λ-4)*x0^(3-λ))/((λ-3)*(λ-4))
    ζx = ζ_x(model,z,Vol,Temp)
    return I*(1-ζx/2)/(1-ζx)^3-9*J*ζx*(ζx+1)/(2*(1-ζx)^3)
end

function x0(model::SAFTVRQMieFamily, z,Vol,Temp,i,j)
    σ           = model.params.sigma[union(i,j)]
    HSd         = (d(model,z,Vol,Temp,i)+d(model,z,Vol,Temp,j))/2
    return σ/HSd
end

function x0eff(model::SAFTVRQMieFamily, z,Vol,Temp,i,j)
    σ_eff       = σeff(model,z,Vol,Temp,i,j)
    HSd         = (d(model,z,Vol,Temp,i)+d(model,z,Vol,Temp,j))/2
    return σ_eff/HSd
end

function a_2(model::SAFTVRQMieFamily, z,Vol,Temp)
    x    = z/sum(z[i] for i in model.components)

    m    = model.params.segment
    m̄    = sum(x[j]*m[j] for j in model.components)
    return m̄/Temp^2*sum(xS(model,z,Vol,Temp,i)*xS(model,z,Vol,Temp,j)*a_2ij(model,z,Vol,Temp,i,j) for i in model.components for j in model.components)
end

function a_2ij(model::SAFTVRQMieFamily, z,Vol,Temp,i,j)
    ϵ           = model.params.epsilon[union(i,j)]
    σ           = model.params.sigma[union(i,j)]
    λR          = model.params.lambdaR[union(i,j)]
    λA          = model.params.lambdaA[union(i,j)]
    Mr          = model.params.MolarMass[union(i,j)]/N_A

    D           = ħ^2/(12*k_B*Temp*Mr*σ^2)

    HSd         = (d(model,z,Vol,Temp,i)+d(model,z,Vol,Temp,j))/2
    x_0eff      = x0eff(model,z,Vol,Temp,i,j)
    x_0         = x0(model,z,Vol,Temp,i,j)
    C           = CMie(model,z,Vol,Temp,i,j)
    ζx          = ζ_x(model,z,Vol,Temp)
    KHS         = (1-ζx)^4/(1+4*ζx+4*ζx^2-4*ζx^3+ζx^4)
    return π*KHS*(1+χ(model,z,Vol,Temp,i,j))*ρs(model,z,Vol,Temp)*ϵ^2*HSd^3*C^2*(x_0^(2*λA)*(a1s(model,z,Vol,Temp,2*λA)+B(model,z,Vol,Temp,2*λA,x_0eff))-
                                                                                 x_0^(λA+λR)*2*(a1s(model,z,Vol,Temp,λA+λR)+B(model,z,Vol,Temp,λA+λR,x_0eff))+
                                                                                 x_0^(2*λR)*(a1s(model,z,Vol,Temp,2*λR)+B(model,z,Vol,Temp,2*λR,x_0eff))+
                                                                                 x_0^(2*λA+2)*2*Q1(λA)*(a1s(model,z,Vol,Temp,2*λA+2)+B(model,z,Vol,Temp,2*λA+2,x_0eff))*D+
                                                                                 x_0^(2*λR+2)*2*Q1(λR)*(a1s(model,z,Vol,Temp,2*λR+2)+B(model,z,Vol,Temp,2*λR+2,x_0eff))*D-
                                                                                 x_0^(λA+λR+2)*2*(Q1(λA)+Q1(λR))*(a1s(model,z,Vol,Temp,λA+λR+2)+B(model,z,Vol,Temp,λA+λR+2,x_0eff))*D+
                                                                                 x_0^(2*λA+4)*Q1(λA)^2*(a1s(model,z,Vol,Temp,2*λA+4)+B(model,z,Vol,Temp,2*λA+4,x_0eff))*D^2+
                                                                                 x_0^(2*λR+4)*Q1(λR)^2*(a1s(model,z,Vol,Temp,2*λR+4)+B(model,z,Vol,Temp,2*λR+4,x_0eff))*D^2-
                                                                                 x_0^(λA+λR+4)*(2*Q1(λA)*Q1(λR))*(a1s(model,z,Vol,Temp,λA+λR+4)+B(model,z,Vol,Temp,λA+λR+4,x_0eff))*D^2+
                                                                                 x_0^(2*λA+4)*2*Q2(λA)*(a1s(model,z,Vol,Temp,2*λA+4)+B(model,z,Vol,Temp,2*λA+4,x_0eff))*D^2+
                                                                                 x_0^(2*λR+4)*2*Q2(λR)*(a1s(model,z,Vol,Temp,2*λR+4)+B(model,z,Vol,Temp,2*λR+4,x_0eff))*D^2-
                                                                                 x_0^(λA+λR+4)*2*(Q2(λA)+Q2(λR))*(a1s(model,z,Vol,Temp,λA+λR+4)+B(model,z,Vol,Temp,λA+λR+4,x_0eff))*D^2+
                                                                                 x_0^(2*λA+6)*2*(Q1(λA)*Q2(λA))*(a1s(model,z,Vol,Temp,2*λA+6)+B(model,z,Vol,Temp,2*λA+6,x_0eff))*D^3+
                                                                                 x_0^(2*λR+6)*2*(Q1(λR)*Q2(λR))*(a1s(model,z,Vol,Temp,2*λR+6)+B(model,z,Vol,Temp,2*λR+6,x_0eff))*D^3-
                                                                                 x_0^(λA+λR+6)*2*(Q1(λR)*Q2(λA)+Q1(λA)*Q2(λR))*(a1s(model,z,Vol,Temp,λA+λR+6)+B(model,z,Vol,Temp,λA+λR+6,x_0eff))*D^3+
                                                                                 x_0^(2*λA+8)*Q2(λA)^2*(a1s(model,z,Vol,Temp,2*λA+8)+B(model,z,Vol,Temp,2*λA+8,x_0eff))*D^4+
                                                                                 x_0^(2*λR+8)*Q2(λR)^2*(a1s(model,z,Vol,Temp,2*λR+8)+B(model,z,Vol,Temp,2*λR+8,x_0eff))*D^4-
                                                                                 x_0^(λA+λR+8)*(2*Q2(λA)*Q2(λR))*(a1s(model,z,Vol,Temp,λA+λR+8)+B(model,z,Vol,Temp,λA+λR+8,x_0eff))*D^4)
end

function χ(model::SAFTVRQMieFamily, z,Vol,Temp,i,j)
    λR          = model.params.lambdaR[union(i,j)]
    λA          = model.params.lambdaA[union(i,j)]
    σ           = model.params.sigma[union(i,j)]
    σ_eff       = σeff(model,z,Vol,Temp,i,j)
    ϵ           = model.params.epsilon[union(i,j)]
    ϵ_eff       = ϵeff(model,z,Vol,Temp,i,j)
    Mr          = model.params.MolarMass[union(i,j)]/N_A

    D           = ħ^2/(12*k_B*Temp*Mr*σ^2)

    ζstar = ζ_star(model, z,Vol,Temp)
    α     = CMie(model,z,Vol,Temp,i,j)*ϵ/ϵ_eff*(((σ/σ_eff)^λA/(λA-3)-(σ/σ_eff)^λR/(λR-3))+
                                              D*((σ/σ_eff)^(2+λA)*Q1(λA)/(λA-1)-(σ/σ_eff)^(2+λR)*Q1(λR)/(λR-1))+
                                            D^2*((σ/σ_eff)^(4+λA)*Q2(λA)/(λA+1)-(σ/σ_eff)^(4+λR)*Q2(λR)/(λR+1)))
    return f(α,1)*ζstar+f(α,2)*ζstar^5+f(α,3)*ζstar^8
end

function ζ_star(model::SAFTVRQMieFamily, z,Vol,Temp)
    return ρs(model,z,Vol,Temp)*π/6*sum(xS(model,z,Vol,Temp,i)*xS(model,z,Vol,Temp,j)*(σeff(model,z,Vol,Temp,i,j))^3 for i in model.components for j in model.components)
end

function a_3(model::SAFTVRQMieFamily, z,Vol,Temp)
    x    = z/sum(z[i] for i in model.components)

    m    = model.params.segment
    m̄    = sum(x[j]*m[j] for j in model.components)
    return m̄/Temp^3*sum(xS(model,z,Vol,Temp,i)*xS(model,z,Vol,Temp,j)*a_3ij(model,z,Vol,Temp,i,j) for i in model.components for j in model.components)
end

function a_3ij(model::SAFTVRQMieFamily, z,Vol,Temp,i,j)
    λR          = model.params.lambdaR[union(i,j)]
    λA          = model.params.lambdaA[union(i,j)]
    σ           = model.params.sigma[union(i,j)]
    σ_eff       = σeff(model,z,Vol,Temp,i,j)
    ϵ           = model.params.epsilon[union(i,j)]
    ϵ_eff       = ϵeff(model,z,Vol,Temp,i,j)
    Mr          = model.params.MolarMass[union(i,j)]/N_A

    D           = ħ^2/(12*k_B*Temp*Mr*σ^2)

    ζstar = ζ_star(model, z,Vol,Temp)
    α     = CMie(model,z,Vol,Temp,i,j)*ϵ/ϵ_eff*(((σ/σ_eff)^λA/(λA-3)-(σ/σ_eff)^λR/(λR-3))+
                                              D*((σ/σ_eff)^(2+λA)*Q1(λA)/(λA-1)-(σ/σ_eff)^(2+λR)*Q1(λR)/(λR-1))+
                                            D^2*((σ/σ_eff)^(4+λA)*Q2(λA)/(λA+1)-(σ/σ_eff)^(4+λR)*Q2(λR)/(λR+1)))
    return -ϵ_eff^3*f(α,4)*ζstar*exp(f(α,5)*ζstar+f(α,6)*ζstar^2)
end
