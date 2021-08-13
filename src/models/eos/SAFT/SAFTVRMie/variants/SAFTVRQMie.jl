struct SAFTVRQMieParam <: EoSParam
    segment::SingleParam{Float64}
    sigma::PairParam{Float64}
    lambda_a::PairParam{Float64}
    lambda_r::PairParam{Float64}
    epsilon::PairParam{Float64}
    Mw::SingleParam{Float64}
end

abstract type SAFTVRQMieModel <: SAFTVRMieModel end
@newmodel SAFTVRQMie SAFTVRQMieModel SAFTVRQMieParam

export SAFTVRQMie
function SAFTVRQMie(components; idealmodel=BasicIdeal, userlocations=String[], ideal_userlocations=String[], verbose=false)
    params = getparams(components, ["SAFT/SAFTVRQMie"]; userlocations=userlocations, verbose=verbose)

    params["Mw"].values .*= 1E-3
    Mw = params["Mw"]
    segment = params["m"]
    params["sigma"].values .*= 1E-10
    sigma = sigma_LorentzBerthelot(params["sigma"])
    epsilon = epsilon_HudsenMcCoubrey(params["epsilon"], sigma)
    lambda_a = lambda_LorentzBerthelot(params["lambda_a"])
    lambda_r = lambda_LorentzBerthelot(params["lambda_r"])

    packagedparams = SAFTVRQMieParam(segment, sigma, lambda_a, lambda_r, epsilon, Mw)
    references = ["todo"]

    model = SAFTVRQMie(packagedparams, idealmodel; ideal_userlocations=ideal_userlocations, references=references, verbose=verbose)
    return model
end


function a_res(model::SAFTVRQMieModel, V, T, z)
    return @f(a_mono)
end

function d(model::SAFTVRQMieModel, V, T, z, i)
    ϵ = model.params.epsilon.diagvalues
    σ = model.params.sigma.diagvalues
    λr = model.params.lambda_r.diagvalues
    λa = model.params.lambda_a.diagvalues
    Mw = model.params.Mw.values

    Di = ħ^2/(12*k_B*T*Mw[i]/N_A*σ[i]^2)
    βi = ϵ[i]/T
    u(x) = @f(C,i,i)*((x^-λr[i]-x^-λa[i])+
              Di*(@f(Q1,λr[i])*x^-(λr[i]+2)-@f(Q1,λa[i])*x^-(λa[i]+2))+
              Di^2*(@f(Q2,λr[i])*x^-(λr[i]+4)-@f(Q2,λa[i])*x^-(λa[i]+4)))
    du(x) = -@f(C,i,i)*((λr[i]*x^-(λr[i]+1)-λa[i]*x^-(λa[i]+1))+
                Di*(@f(Q1,λr[i])*(λr[i]+2)*x^-(λr[i]+3)-@f(Q1,λa[i])*(λa[i]+2)*x^-(λa[i]+3))+
                Di^2*(@f(Q2,λr[i])*(λr[i]+4)*x^-(λr[i]+5)-@f(Q2,λa[i])*(λa[i]+4)*x^-(λa[i]+5)))
    d2u(x) = @f(C,i,i)*((λr[i]*(λr[i]+1)*x^-(λr[i]+2)-λa[i]*(λa[i]+1)*x^-(λa[i]+2))+
                Di*(@f(Q1,λr[i])*(λr[i]+2)*(λr[i]+3)*x^-(λr[i]+4)-@f(Q1,λa[i])*(λa[i]+2)*(λa[i]+3)*x^-(λa[i]+4))+
                Di^2*(@f(Q2,λr[i])*(λr[i]+4)*(λr[i]+5)*x^-(λr[i]+6)-@f(Q2,λa[i])*(λa[i]+4)*(λa[i]+5)*x^-(λa[i]+6)))

    f(x) = exp(-βi*u(x))
    g(x) = -βi*exp(-βi*u(x))*du(x)
    h(x) = βi*exp(-βi*u(x))*(βi*du(x)^2-d2u(x))
    x_min = @f(Halley,f,g,h,1.)
    σ_effi = @f(σeff,i,i)/σ[i]
    x = SAFTVRQMieconsts.x
    w = SAFTVRQMieconsts.w
    return σ[i]*(σ_effi-(σ_effi-x_min)/2*∑(w[i]*f((σ_effi-x_min)/2*x[i]+(σ_effi+x_min)/2) for i in 1:21))
end

function σeff(model::SAFTVRQMieModel, V, T, z, i, j)
    σ = model.params.sigma.values
    λr = model.params.lambda_r.values
    λa = model.params.lambda_a.values
    Mwij = (model.params.Mw.values[i] + model.params.Mw.values[j])/2 # check
    Dij = ħ^2/(12*k_B*T*Mwij/N_A*σ[i,j]^2)

    f(x) = Dij^2*@f(Q2,λr[i,j])+Dij*@f(Q1,λr[i,j])*x^2+x^4-x^(λr[i,j]-λa[i,j])*(Dij^2*@f(Q2,λa[i,j])+Dij*@f(Q1,λa[i,j])*x^2+x^4)
    g(x) = 2*Dij*@f(Q1,λr[i,j])*x+4x^3-x^(λr[i,j]-λa[i,j]-1)*((λr[i,j]-λa[i,j])*Dij^2*@f(Q2,λa[i,j])+(λr[i,j]-λa[i,j]+2)*Dij*@f(Q1,λa[i,j])*x^2+(λr[i,j]-λa[i,j]+4)*x^4)
    h(x) = 2*Dij*@f(Q1,λr[i,j])+12x^2-x^(λr[i,j]-λa[i,j]-2)*((λr[i,j]-λa[i,j]-1)*(λr[i,j]-λa[i,j])*Dij^2*@f(Q2,λa[i,j])+((λr[i,j]-λa[i,j]-1)*(λr[i,j]-λa[i,j]+2)+2*(λr[i,j]-λa[i,j]+2))*Dij*@f(Q1,λa[i,j])*x^2+((λr[i,j]-λa[i,j]-1)*(λr[i,j]-λa[i,j]+4)+4*(λr[i,j]-λa[i,j]+4))*x^4)

    return σ[i,j]*@f(Halley,f,g,h,1.)
end

function ϵeff(model::SAFTVRQMieModel, V, T, z, i, j)
    ϵ = model.params.epsilon.values
    σ = model.params.sigma.values
    λr = model.params.lambda_r.values
    λa = model.params.lambda_a.values
    Mwij = (model.params.Mw.values[i] + model.params.Mw.values[j])/2 # check
    Dij = ħ^2/(12*k_B*T*Mwij/N_A*σ[i,j]^2)

    u(x) = @f(C,i,j)*((x^-λr[i,j]-x^-λa[i,j])+
            Dij*(@f(Q1,λr[i,j])*x^-(λr[i,j]+2)-@f(Q1,λa[i,j])*x^-(λa[i,j]+2))+
            Dij^2*(@f(Q2,λr[i,j])*x^-(λr[i,j]+4)-@f(Q2,λa[i,j])*x^-(λa[i,j]+4)))
    f(x) = (λr[i,j]+4)*Dij^2*@f(Q2,λr[i,j])+(λr[i,j]+2)*Dij*@f(Q1,λr[i,j])*x^2+λr[i,j]*x^4-x^(λr[i,j]-λa[i,j])*((λa[i,j]+4)*Dij^2*@f(Q2,λa[i,j])+(λa[i,j]+2)*Dij*@f(Q1,λa[i,j])*x^2+λa[i,j]*x^4)
    g(x) = 2*(λr[i,j]+2)*Dij*@f(Q1,λr[i,j])*x+4λr[i,j]*x^3-x^(λr[i,j]-λa[i,j]-1)*((λa[i,j]+4)*(λr[i,j]-λa[i,j])*Dij^2*@f(Q2,λa[i,j])+(λa[i,j]+2)*(λr[i,j]-λa[i,j]+2)*Dij*@f(Q1,λa[i,j])*x^2+λa[i,j]*(λr[i,j]-λa[i,j]+4)*x^4)
    h(x) = 2*(λr[i,j]+2)*Dij*@f(Q1,λr[i,j])+12λr[i,j]*x^2-x^(λr[i,j]-λa[i,j]-2)*((λa[i,j]+4)*(λr[i,j]-λa[i,j]-1)*(λr[i,j]-λa[i,j])*Dij^2*@f(Q2,λa[i,j])+(λa[i,j]+2)*((λr[i,j]-λa[i,j]-1)*(λr[i,j]-λa[i,j]+2)+2*(λr[i,j]-λa[i,j]+2))*Dij*@f(Q1,λa[i,j])*x^2+λa[i,j]*((λr[i,j]-λa[i,j]-1)*(λr[i,j]-λa[i,j]+4)+4*(λr[i,j]-λa[i,j]+4))*x^4)

    σ_min = @f(Halley,f,g,h,(λr[i,j]/λa[i,j])^(1/(λr[i,j]-λa[i,j])))
    return -ϵ[i,j]*u(σ_min)
end

function Halley(model::SAFTVRQMieModel, V, T, z, f_, g_, h_, x0)
    tolx = 1.
    tolf = 1.
    f0 = f_(x0)
    if f0 < 1e-16
        return x0
    else
        while tolf > 1e-16 && tolx > 1e-15
            d = f_(x0)/g_(x0)*(1-f_(x0)*h_(x0)/(2*g_(x0)^2))^-1
            x0=x0-d
            f0 = f_(x0)
            tolf = abs(f0)
            tolx = abs(d)
        end
        return x0
    end
end

function ζ_x(model::SAFTVRQMieModel, V, T, z)
    comps = @comps
    return π/6*@f(ρ_S)*∑(@f(x_S,i)*@f(x_S,j)*(@f(d,i,j))^3 for i in comps for j in comps)
end

function a_1(model::SAFTVRQMieModel, V, T, z, i, j)
    ϵ = model.params.epsilon.values
    σ = model.params.sigma.values
    λr = model.params.lambda_r.values
    λa = model.params.lambda_a.values
    Mwij = (model.params.Mw.values[i] + model.params.Mw.values[j])/2 # check
    Dij = ħ^2/(12*k_B*T*Mwij/N_A*σ[i,j]^2)
    x_0ij = @f(x_0,i,j)
    x_0effij = @f(x_0eff,i,j)
    return 2*π*ϵ[i,j]*@f(d,i,j)^3*@f(C,i,j)*@f(ρ_S)*
    ( (x_0ij^λa[i,j]*(@f(aS_1,λa[i,j])+@f(B,λa[i,j],x_0effij))-
        x_0ij^λr[i,j]*(@f(aS_1,λr[i,j])+@f(B,λr[i,j],x_0effij)))+
        (x_0ij^(λa[i,j]+2)*@f(Q1,λa[i,j])*(@f(aS_1,λa[i,j]+2)+@f(B,λa[i,j]+2,x_0effij))-
         x_0ij^(λr[i,j]+2)*@f(Q1,λr[i,j])*(@f(aS_1,λr[i,j]+2)+@f(B,λr[i,j]+2,x_0effij)))*Dij+
        (x_0ij^(λa[i,j]+4)*@f(Q2,λa[i,j])*(@f(aS_1,λa[i,j]+4)+@f(B,λa[i,j]+4,x_0effij))-
         x_0ij^(λr[i,j]+4)*@f(Q2,λr[i,j])*(@f(aS_1,λr[i,j]+4)+@f(B,λr[i,j]+4,x_0effij)))*Dij^2 )
end

function Q1(model::SAFTVRQMieModel, V, T, z, λ)
    return λ*(λ-1)
end

function Q2(model::SAFTVRQMieModel, V, T, z, λ)
    return 1/2*(λ+2)*(λ+1)*λ*(λ-1)
end

function x_0eff(model::SAFTVRQMieModel, V, T, z, i, j)
    return @f(σeff,i,j)/@f(d,i,j)
end

function a_2(model::SAFTVRQMieModel, V, T, z, i, j)
    ϵ = model.params.epsilon.values
    σ = model.params.sigma.values
    λr = model.params.lambda_r.values[i,j]
    λa = model.params.lambda_a.values[i,j]
    Mwij = (model.params.Mw.values[i] + model.params.Mw.values[j])/2 # check
    Dij = ħ^2/(12*k_B*T*Mwij/N_A*σ[i,j]^2)

    x_0ij = @f(x_0,i,j)
    x_0effij = @f(x_0eff,i,j)
    return π*@f(KHS)*(1+@f(χ,i,j))*@f(ρ_S)*ϵ[i,j]^2*@f(d,i,j)^3*@f(C,i,j)^2*(x_0ij^(2*λa)*(@f(aS_1,2*λa)+@f(B,2*λa,x_0effij))-
    x_0ij^(λa+λr)*2*(@f(aS_1,λa+λr)+@f(B,λa+λr,x_0effij))+
    x_0ij^(2*λr)*(@f(aS_1,2*λr)+@f(B,2*λr,x_0effij))+
    x_0ij^(2*λa+2)*2*@f(Q1,λa)*(@f(aS_1,2*λa+2)+@f(B,2*λa+2,x_0effij))*Dij+
    x_0ij^(2*λr+2)*2*@f(Q1,λr)*(@f(aS_1,2*λr+2)+@f(B,2*λr+2,x_0effij))*Dij-
    x_0ij^(λa+λr+2)*2*(@f(Q1,λa)+@f(Q1,λr))*(@f(aS_1,λa+λr+2)+@f(B,λa+λr+2,x_0effij))*Dij+
    x_0ij^(2*λa+4)*@f(Q1,λa)^2*(@f(aS_1,2*λa+4)+@f(B,2*λa+4,x_0effij))*Dij^2+
    x_0ij^(2*λr+4)*@f(Q1,λr)^2*(@f(aS_1,2*λr+4)+@f(B,2*λr+4,x_0effij))*Dij^2-
    x_0ij^(λa+λr+4)*(2*@f(Q1,λa)*@f(Q1,λr))*(@f(aS_1,λa+λr+4)+@f(B,λa+λr+4,x_0effij))*Dij^2+
    x_0ij^(2*λa+4)*2*@f(Q2,λa)*(@f(aS_1,2*λa+4)+@f(B,2*λa+4,x_0effij))*Dij^2+
    x_0ij^(2*λr+4)*2*@f(Q2,λr)*(@f(aS_1,2*λr+4)+@f(B,2*λr+4,x_0effij))*Dij^2-
    x_0ij^(λa+λr+4)*2*(@f(Q2,λa)+@f(Q2,λr))*(@f(aS_1,λa+λr+4)+@f(B,λa+λr+4,x_0effij))*Dij^2+
    x_0ij^(2*λa+6)*2*(@f(Q1,λa)*@f(Q2,λa))*(@f(aS_1,2*λa+6)+@f(B,2*λa+6,x_0effij))*Dij^3+
    x_0ij^(2*λr+6)*2*(@f(Q1,λr)*@f(Q2,λr))*(@f(aS_1,2*λr+6)+@f(B,2*λr+6,x_0effij))*Dij^3-
    x_0ij^(λa+λr+6)*2*(@f(Q1,λr)*@f(Q2,λa)+@f(Q1,λa)*@f(Q2,λr))*(@f(aS_1,λa+λr+6)+@f(B,λa+λr+6,x_0effij))*Dij^3+
    x_0ij^(2*λa+8)*@f(Q2,λa)^2*(@f(aS_1,2*λa+8)+@f(B,2*λa+8,x_0effij))*Dij^4+
    x_0ij^(2*λr+8)*@f(Q2,λr)^2*(@f(aS_1,2*λr+8)+@f(B,2*λr+8,x_0effij))*Dij^4-
    x_0ij^(λa+λr+8)*(2*@f(Q2,λa)*@f(Q2,λr))*(@f(aS_1,λa+λr+8)+@f(B,λa+λr+8,x_0effij))*Dij^4)
end

function χ(model::SAFTVRQMieModel, V, T, z,i,j)
    λr = model.params.lambda_r.values
    λa = model.params.lambda_a.values
    σ = model.params.sigma.values
    ϵ = model.params.epsilon.values
    σeffij = @f(σeff,i,j)
    ϵeffij = @f(ϵeff,i,j)
    Mwij = (model.params.Mw.values[i] + model.params.Mw.values[j])/2 # check
    Dij = ħ^2/(12*k_B*T*Mwij/N_A*σ[i,j]^2)

    ζst_ = ζst(model, V, T, z)
    α = @f(C,i,j)*ϵ[i,j]/ϵeffij*
    (((σ[i,j]/σeffij)^λa[i,j]/(λa[i,j]-3)-(σ[i,j]/σeffij)^λr[i,j]/(λr[i,j]-3))+
        Dij*((σ[i,j]/σeffij)^(2+λa[i,j])*@f(Q1,λa[i,j])/(λa[i,j]-1)-
             (σ[i,j]/σeffij)^(2+λr[i,j])*@f(Q1,λr[i,j])/(λr[i,j]-1))+
        Dij^2*((σ[i,j]/σeffij)^(4+λa[i,j])*@f(Q2,λa[i,j])/(λa[i,j]+1)-
               (σ[i,j]/σeffij)^(4+λr[i,j])*@f(Q2,λr[i,j])/(λr[i,j]+1)))
    return @f(f,α,1)*ζst_+@f(f,α,2)*ζst_^5+@f(f,α,3)*ζst_^8
end

function ζst(model::SAFTVRQMieModel, V, T, z)
    comps = @comps
    return @f(ρ_S)*π/6*∑(@f(x_S,i)*@f(x_S,j)*(@f(σeff,i,j))^3 for i in comps for j in comps)
end

function a_3(model::SAFTVRQMieModel, V, T, z, i, j)
    λr = model.params.lambda_r.values
    λa = model.params.lambda_a.values
    σ = model.params.sigma.values
    ϵ = model.params.epsilon.values
    σeffij = @f(σeff,i,j)
    ϵeffij = @f(ϵeff,i,j)
    Mwij = (model.params.Mw.values[i] + model.params.Mw.values[j])/2 # check
    Dij = ħ^2/(12*k_B*T*Mwij/N_A*σ[i,j]^2)

    ζst_ = ζst(model, V, T, z)
    α = @f(C,i,j)*ϵ[i,j]/ϵeffij*
    (((σ[i,j]/σeffij)^λa[i,j]/(λa[i,j]-3)-(σ[i,j]/σeffij)^λr[i,j]/(λr[i,j]-3))+
        Dij*((σ[i,j]/σeffij)^(2+λa[i,j])*@f(Q1,λa[i,j])/(λa[i,j]-1)-
             (σ[i,j]/σeffij)^(2+λr[i,j])*@f(Q1,λr[i,j])/(λr[i,j]-1))+
        Dij^2*((σ[i,j]/σeffij)^(4+λa[i,j])*@f(Q2,λa[i,j])/(λa[i,j]+1)-
               (σ[i,j]/σeffij)^(4+λr[i,j])*@f(Q2,λr[i,j])/(λr[i,j]+1)))
    return -ϵeffij^3*@f(f,α,4)*ζst_*exp(@f(f,α,5)*ζst_+@f(f,α,6)*ζst_^2)
end

const SAFTVRQMieconsts = (
    x = [-0.9956571630258080807355,-0.973906528517171720078,-0.9301574913557082260012,-0.8650633666889845107321,-0.7808177265864168970637,-0.6794095682990244062343,-0.562757134668604683339,
        -0.4333953941292471907993,-0.294392862701460198131,-0.1488743389816312108848,0,0.1488743389816312108848,0.2943928627014601981311,0.4333953941292471907993,
        0.562757134668604683339,0.6794095682990244062343,0.7808177265864168970637,0.865063366688984510732,0.9301574913557082260012,0.973906528517171720078,0.9956571630258080807355],
    w = [0.0116946388673718742781,0.0325581623079647274788,0.0547558965743519960314,0.075039674810919952767,0.093125454583697605535,0.1093871588022976418992,0.123491976262065851078,
        0.134709217311473325928,0.142775938577060080797,0.1477391049013384913748,0.149445554002916905665,0.1477391049013384913748,0.1427759385770600807971,0.134709217311473325928,
        0.123491976262065851078,0.109387158802297641899,0.093125454583697605535,0.075039674810919952767,0.05475589657435199603138,0.032558162307964727479,0.0116946388673718742781],
)
