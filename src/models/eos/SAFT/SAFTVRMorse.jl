const SAFTVRMorseconsts = (
    u = [-0.973906528517172,-0.865063366688985,-0.679409568299024,-0.433395394129247,-0.148874338981631,0.148874338981631,0.433395394129247,0.679409568299024,0.865063366688985,0.973906528517172]
    ,w = [0.0666713443086881,0.149451349150581,0.219086362515982,0.269266719309996,0.295524224714753,0.295524224714753,0.269266719309996,0.219086362515982,0.149451349150581,0.0666713443086881]

    ,τ1 = [-0.005929462731  0.028821220185  0.823222580046 -1.514002220875 -25.566034300040;
            0.004297313190 -0.004085322342 -0.655058954616 -1.434307204889  37.285171618599;
           -0.001492689050  0.023497433501 -0.235415060297  3.380763792236 -18.166122057226;
           -0.000498507895  0.016751329671 -0.161840634134  0.449938407721  -0.293219371958;
            0.000080667308 -0.000405514821 -0.028927870642  0.452624568510  -2.226513843436]
    ,τ2 = [ 0.023134461 -0.57548818   6.490532003 -35.1940007   69.11113346;
            0.007072932 -0.216013354  1.890413397 -6.002948117  3.222116026;
            0.007901316 -0.179042936  1.404582995 -3.823782235  2.858870428;
            0.002074109 -0.055995333  0.622738789 -3.366070495  7.445109302;
            0.000958321 -0.022654207  0.206722405 -0.917821835  2.041301543]
    ,τchain = [ -0.368126173 -0.005107786  -0.048071831 -0.022541511;
                  4.253900887  0.653256069   0.482204788  0.374467426;
                -24.852056185 -5.552055928  -2.074631541 -2.961697758;
                 38.899991521  11.810423426  4.032172895 10.557991805]

    ,ϕ = [18.35352712    7552.19773298   3902503.69092418;
          -4.53981344    45122.82036360  11294513.73984250;
          -113.72796238  75074.64229911 -507562.17665894;
          -118.08065115  38471.46257363  8832404.17580299;
          -58.26532521   0.27323257     -4956.05044240;
          -9.06145956   -15.12866683     3061986.36256341;
          -41.46384272  -23.40145163     1310.11514429;
          -14.59649156  -56.38187762     413.27679449]

    ,c  = [0.0756425183020431	-0.128667137050961	 0.128350632316055	-0.0725321780970292	   0.0257782547511452  -0.00601170055221687	  0.000933363147191978  -9.55607377143667e-05  6.19576039900837e-06 -2.30466608213628e-07 3.74605718435540e-09
          0.134228218276565	    -0.182682168504886 	 0.0771662412959262	-0.000717458641164565 -0.00872427344283170	0.00297971836051287	 -0.000484863997651451	 4.35262491516424e-05 -2.07789181640066e-06	4.13749349344802e-08 0
         -0.565116428942893	     1.00930692226792   -0.660166945915607	 0.214492212294301	  -0.0388462990166792	0.00406016982985030	 -0.000239515566373142	 7.25488368831468e-06 -8.58904640281928e-08	0	                 0
         -0.387336382687019	    -0.211614570109503	 0.450442894490509	-0.176931752538907	   0.0317171522104923  -0.00291368915845693	  0.000130193710011706  -2.14505500786531e-06  0	                0	                 0
          2.13713180911797	    -2.02798460133021 	 0.336709255682693	 0.00118106507393722  -0.00600058423301506	0.000626343952584415 -2.03636395699819e-05	 0	                   0	                0	                 0
         -0.300527494795524	     2.89920714512243   -0.567134839686498	 0.0518085125423494	  -0.00239326776760414	4.15107362643844e-05  0	                     0	                   0	                0                    0
         -6.21028065719194	    -1.92883360342573	 0.284109761066570	-0.0157606767372364	   0.000368599073256615	0 	                  0	                     0	                   0	                0	                 0
          11.6083532818029	     0.742215544511197  -0.0823976531246117	 0.00186167650098254   0	                0	                  0	                     0	                   0	                0	                 0
         -10.2632535542427	    -0.125035689035085	 0.0114299144831867	 0	                   0	                0	                  0	                     0	                   0	                0	                 0
          4.65297446837297	    -0.00192518067137033 0	                 0	                   0	                0	                  0	                     0	                   0	                0	                 0
         -0.867296219639940	     0	                 0	                 0	                   0	                0	                  0	                     0	                   0	                0	                 0]
)

function a_res(model::SAFTVRMorseFamily, z,Vol,Temp)
    return a_mono(model,z,Vol,Temp)+a_chain(model,z,Vol,Temp)+a_assoc(model,z,Vol,Temp)
end

function a_mono(model::SAFTVRMorseFamily, z,Vol,Temp)
    return a_hs(model,z,Vol,Temp)+a_disp(model,z,Vol,Temp)
end
function a_disp(model::SAFTVRMorseFamily, z,Vol,Temp)
    return a_1(model,z,Vol,Temp)+a_2(model,z,Vol,Temp)
end

function a_hs(model::SAFTVRMorseFamily, z,Vol,Temp)
    ζ0   = ζn(model, z,Vol,Temp, 0)
    ζ1   = ζn(model, z,Vol,Temp, 1)
    ζ2   = ζn(model, z,Vol,Temp, 2)
    ζ3   = ζn(model, z,Vol,Temp, 3)
    NParticles = N_A*sum(z[i] for i in model.components)
    return 6*Vol/π/NParticles*(3ζ1*ζ2/(1-ζ3) + ζ2^3/(ζ3*(1-ζ3)^2) + (ζ2^3/ζ3^2-ζ0)*log(1-ζ3))
end

function ζn(model::SAFTVRMorseFamily, z,Vol,Temp,n)
    return π/6*ρs(model,z,Vol,Temp)*sum(xS(model,z,Vol,Temp,i)*d(model,z,Vol,Temp,i)^n for i in model.components)
end

function ρs(model::SAFTVRMorseFamily, z,Vol,Temp)
    NParticles = N_A*sum(z[i] for i in model.components)
    x   = z/sum(z[i] for i in model.components)
    m   = model.params.segment
    m̄   = sum(x[i]*m[i] for i in model.components)
    return NParticles/Vol*m̄
end

function xS(model::SAFTVRMorseFamily, z,Vol,Temp,i)
    x   = z/sum(z[i] for i in model.components)
    m   = model.params.segment
    m̄   = sum(x[j]*m[j] for j in model.components)
    return x[i]*m[i]/m̄
end

function d(model::SAFTVRMorseFamily, z,Vol,Temp,i)
    ϵ           = model.params.epsilon[union(i,i)]
    σ           = model.params.sigma[union(i,i)]
    λ           = model.params.lambda[union(i,i)]
    u           = SAFTVRMorseconsts.u
    w           = SAFTVRMorseconsts.w
    T̄           = Temp/ϵ
    ϕ(x)        = (exp(-2λ*(x-2^(1/6)))-2*exp(-λ*(x-2^(1/6))))
    return σ*(1-sum(w[i]*exp(-ϕ(u[i]/2+1/2)/T̄) for i in 1:10)/2)
end

function a_1(model::SAFTVRMorseFamily, z,Vol,Temp)
    x    = z/sum(z[i] for i in model.components)

    m    = model.params.segment
    m̄    = sum(x[j]*m[j] for j in model.components)
    return 2π*N_A*sum(z)/Vol*m̄^2/Temp*sum(xS(model,z,Vol,Temp,i)*xS(model,z,Vol,Temp,j)*a_1ij(model,z,Vol,Temp,i,j) for i in model.components for j in model.components)
end

function a_1ij(model::SAFTVRMorseFamily, z,Vol,Temp,i,j)
    ϵ           = model.params.epsilon[union(i,j)]
    λ           = model.params.lambda[union(i,j)]

    HSd         = (d(model,z,Vol,Temp,i)+d(model,z,Vol,Temp,j))/2
    η           = ζn(model,z,Vol,Temp,3)
    return ϵ*HSd^3*I1(model,z,Vol,Temp,λ,η)
end

function I1(model::SAFTVRMorseFamily,z,Vol,Temp,λ,η)
    return sum(c1(model,z,Vol,Temp,λ,n)*η^(5-n) for n in 1:5)
end

function c1(model::SAFTVRMorseFamily,z,Vol,Temp,λ,n)
    τ1 = SAFTVRMorseconsts.τ1
    return sum(τ1[n,i]*λ^(5-i) for i in 1:5)
end

function a_2(model::SAFTVRMorseFamily, z,Vol,Temp)
    x    = z/sum(z[i] for i in model.components)

    m    = model.params.segment
    m̄    = sum(x[j]*m[j] for j in model.components)
    ζx          = ζn(model,z,Vol,Temp,3)
    KHS         = (1-ζx)^4/(1+4*ζx+4*ζx^2-4*ζx^3+ζx^4)
    return π*N_A*sum(z)/Vol*m̄^2/Temp^2*sum(xS(model,z,Vol,Temp,i)*xS(model,z,Vol,Temp,j)*a_2ij(model,z,Vol,Temp,i,j) for i in model.components for j in model.components)
end

function a_2ij(model::SAFTVRMorseFamily, z,Vol,Temp,i,j)
    ϵ           = model.params.epsilon[union(i,j)]
    λ           = model.params.lambda[union(i,j)]
    σ           = model.params.sigma[union(i,j)]

    HSd         = (d(model,z,Vol,Temp,i)+d(model,z,Vol,Temp,j))/2
    x0          = σ/HSd
    η           = ζn(model,z,Vol,Temp,3)
    return -ϵ^2*HSd^3*(1+χ(model,z,Vol,Temp,λ,σ,η,x0))*I2(model,z,Vol,Temp,λ,η)
end

function I2(model::SAFTVRMorseFamily,z,Vol,Temp,λ,η)
    return sum(c2(model,z,Vol,Temp,λ,n)*η^(5-n) for n in 1:5)
end

function c2(model::SAFTVRMorseFamily,z,Vol,Temp,λ,n)
    τ2 = SAFTVRMorseconsts.τ2
    return sum(τ2[n,i]*λ^(5-i) for i in 1:5)
end

function χ(model::SAFTVRMorseFamily, z,Vol,Temp,λ,σ,η,x0)
    re = 2^(1/6)
    α = ((2λ*(λ+1)+1)*exp(2λ*(re-1)))/(4λ^3)-((2λ*(λ+2)+2)*exp(λ*(re-1)))/(λ^3)
    return f(model,α,1)*(x0^3*η)+f(model,α,2)*(x0^3*η)^5+f(model,α,3)*(x0^3*η)^8
end

function f(model::SAFTVRMorse,α,m)
    ϕ  = SAFTVRMorseconsts.ϕ
    return sum(ϕ[i+1,m]*α^i for i in 0:4)/(1+sum(ϕ[i+1,m]*α^(i-3) for i in 5:7))
end

function a_chain(model::SAFTVRMorseFamily, z,Vol,Temp)
    x       = z/sum(z[i] for i in model.components)
    m       = model.params.segment
    return -sum(x[i]*(log(g_Morse(model,z,Vol,Temp,i))*(m[i]-1)) for i in model.components)
end
#
function g_Morse(model::SAFTVRMorseFamily,z,Vol,Temp,i)
    ϵ    = model.params.epsilon[union(i,i)]
    gHS  = g_HS(model,z,Vol,Temp,i)
    g1   = g_1(model,z,Vol,Temp,i)
    gMie = gHS*exp(ϵ/Temp*g1/gHS);
    return gMorse
end
#
function g_HS(model::SAFTVRMorseFamily,z,Vol,Temp,i)
    η   = ζn(model,z,Vol,Temp,3)
    x_0 = model.params.sigma[i]/d(model,z,Vol,Temp,i)
    k0  = -log(1-η)+(42η-39η^2+9η^3-2η^4)/(6*(1-η)^3)
    k1  = (η^4+6*η^2-12*η)/(2*(1-η)^3)
    k2  = -3*η^2/(8*(1-η)^2)
    k3  = (-η^4+3*η^2+3*η)/(6*(1-η)^3)
    gHS = exp(k0+x_0*k1+x_0^2*k2+x_0^3*k3)
    return gHS
end

function g_1(model::SAFTVRMorseFamily,z,Vol,Temp,i)
    λ           = model.params.lambda[i]
    η           = ζn(model,z,Vol,Temp,3)
    g1 = da1(model,z,Vol,Temp,i)-Ichain(model,z,Vol,Temp,λ,η)
    return g1
end

function da1(model::SAFTVRMorseFamily,z,Vol,Temp,i)
    x           = z/sum(z[i] for i in model.components)

    ϵ           = model.params.epsilon[j]
    λ           = model.params.lambda[i]
    m̄           = sum(x[i]*model.params.segment[i] for i in model.components)

    HSd         = d(model,z,Vol,Temp,i)
    η           = ζn(model,z,Vol,Temp,3)
    return (1/Temp)/4*HSd^3*(ρs(model,z,Vol,Temp)/η*I1(model,z,Vol,Temp,λ,η)+ρs(model,z,Vol,Temp)*dI1(model,z,Vol,Temp,λ,η))
end

function dI1(model::SAFTVRMorseFamily,z,Vol,Temp,λ,η)
    return sum((5-n)*c1(model,z,Vol,Temp,λ,n)*η^(4-n) for n in 1:4)
end

function Ichain(model::SAFTVRMorseFamily,z,Vol,Temp,λ,η)
    return sum(c_chain(model,z,Vol,Temp,λ,n)*η^(4-n) for n in 1:4)
end

function c_chain(model,z,Vol,Temp,λ,n)
    τchain = SAFTVRMorseconsts.τchain
    return sum(τchain[n,i]*λ^(4-i) for i in 1:4)
end


function a_assoc(model::SAFTVRMorseFamily, z, v, T)
    x = z/sum(z[i] for i in model.components)
    n_sites = model.params.n_sites
    X_iA = X_assoc(model,z,v,T)
    return sum(x[i]*sum(n_sites[i][a]*(log(X_iA[i,a])+(1-X_iA[i,a])/2) for a in keys(model.params.n_sites[i])) for i in model.components)
end

function X_assoc(model::SAFTVRMorseFamily, z, v, T)
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

function Δ(model::SAFTVRMorseFamily, z, v, T, i, j, a, b)
    ρR = ρs(model,z,v,T)*σx3(model,z,v,T)
    TR = T/model.params.epsilon[union(i,j)]
    σ  = model.params.sigma[union(i,j)]

    c  = SAFTVRMorseconsts.c
    I = sum(sum(c[n+1,m+1]*ρR^n*TR^m for m in 0:(10-n)) for n in 0:10)
    ϵ_assoc = model.params.epsilon_assoc[Set([(i,a),(j,b)])]
    F = (exp(ϵ_assoc/T)-1)
    K = model.params.bond_vol[Set([(i,a),(j,b)])]*σ^3
    return F*K*I
end

function σx3(model::SAFTVRMorseFamily, z, v, T)
    σ = model.params.sigma
    return sum(sum(xS(model,z,v,T,i)*xS(model,z,v,T,j)*σ[union(i,j)]^3 for j in model.components) for i in model.components)
end
