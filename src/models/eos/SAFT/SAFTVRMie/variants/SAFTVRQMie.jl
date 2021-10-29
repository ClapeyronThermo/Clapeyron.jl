struct SAFTVRQMieParam <: EoSParam
    segment::SingleParam{Float64}
    sigma::PairParam{Float64}
    lambda_a::PairParam{Float64}
    lambda_r::PairParam{Float64}
    epsilon::PairParam{Float64}
    Mw::PairParam{Float64}
end

abstract type SAFTVRQMieModel <: SAFTVRMieModel end
@newmodel SAFTVRQMie SAFTVRQMieModel SAFTVRQMieParam

export SAFTVRQMie
function SAFTVRQMie(components; idealmodel=BasicIdeal, userlocations=String[], ideal_userlocations=String[], verbose=false)
    params = getparams(components, ["SAFT/SAFTVRQMie"]; userlocations=userlocations, verbose=verbose)

    params["Mw"].values .*= 1E-3
    Mw = sigma_LorentzBerthelot(params["Mw"])
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

function a_mono(model::SAFTVRQMieModel, V, T, z,_data = @f(data))
    _,_,vrdata = _data
    return @f(a_hs,vrdata) + @f(a_disp,_data)
end

function a_res(model::SAFTVRQMieModel, V, T, z)
    @f(a_mono)
end

function data(model::SAFTVRQMieModel, V, T, z)
    _σeff = @f(σeff)
    _ϵff = @f(ϵeff)
    _d = @f(d,_σeff)
    ζi = @f(ζ0123,_d)
    _ζst = @f(ζst,_σeff)
    _ζ_X,_ = @f(ζ_X_σ3,_d)
    _ρ_S = @f(ρ_S)
    σ3x = _ζst/(_ρ_S*π/6)
    vrdata = (_d,_ρ_S,ζi,_ζ_X,_ζst,σ3x)
    return (_σeff,_ϵff,vrdata)
end

function Q1(model::SAFTVRQMieModel, V, T, z, λ)
    return λ*(λ-1)
end

function Q2(model::SAFTVRQMieModel, V, T, z, λ)
    return 1/2*(λ+2)*(λ+1)*λ*(λ-1)
end

function d(model::SAFTVRQMieModel, V, T, z, i,_data)
    ϵ,σ,λr,λa,Mw,σeff = _data

    Di = ħ^2/(12*k_B*T*Mw/N_A*σ^2)
    βi = ϵ/T

    function udud2u(x)
        _C = @f(Cλ,λa,λr)
        Q1λr = @f(Q1,λr)
        Q1λa = @f(Q1,λa)
        Q2λr = @f(Q2,λr)
        Q2λa = @f(Q2,λa)
        _u= _C*((x^-λr-x^-λa)+
            Di*(Q1λr*x^-(λr+2)-Q1λa*x^-(λa+2))+
            Di^2*(Q2λr*x^-(λr+4)-Q2λa*x^-(λa+4)))
        _du = -_C*((λr*x^-(λr+1)-λa*x^-(λa+1))+
            Di*(Q1λr*(λr+2)*x^-(λr+3)-Q1λa*(λa+2)*x^-(λa+3))+
            Di^2*(Q2λr*(λr+4)*x^-(λr+5)-Q2λa*(λa+4)*x^-(λa+5)))
        _d2u = _C*((λr*(λr+1)*x^-(λr+2)-λa*(λa+1)*x^-(λa+2))+
            Di*(Q1λr*(λr+2)*(λr+3)*x^-(λr+4)-Q1λa*(λa+2)*(λa+3)*x^-(λa+4))+
            Di^2*(Q2λr*(λr+4)*(λr+5)*x^-(λr+6)-Q2λa*(λa+4)*(λa+5)*x^-(λa+6)))
        
        return _u,_du,_d2u
    end

    function f(x)
        _C = @f(Cλ,λa,λr)
        Q1λr = @f(Q1,λr)
        Q1λa = @f(Q1,λa)
        Q2λr = @f(Q2,λr)
        Q2λa = @f(Q2,λa)
        _u= _C*((x^-λr-x^-λa)+
            Di*(Q1λr*x^-(λr+2)-Q1λa*x^-(λa+2))+
            Di^2*(Q2λr*x^-(λr+4)-Q2λa*x^-(λa+4)))
        return exp(-βi*_u)
    end

    function fgh(x)
        _u,_du,_d2u = udud2u(x)
        _f = exp(-βi*_u)
        _g = -βi*_f*_du
        _h = βi*_f*(βi*_du^2-_d2u)
        return _f,_g,_h
    end
    x_min = Solvers.halley(fgh,one(T))
    σ_effi = σeff/σ
    x = SAFTVRQMieconsts.x
    w = SAFTVRQMieconsts.w
    return σ*(σ_effi-(σ_effi-x_min)/2*sum(w[i]*f((σ_effi-x_min)/2*x[i]+(σ_effi+x_min)/2) for i in 1:21))
end

function d(model::SAFTVRQMieModel, V, T, z,_σeff = @f(σeff))
    _ϵ = model.params.epsilon.values
    _σ = model.params.sigma.values
    _λr = model.params.lambda_r.values
    _λa = model.params.lambda_a.values
    _Mwij = model.params.Mw.values
    
    _d = zeros(typeof(T),length(z))
    for i ∈ 1:length(_d)
        _data = (_ϵ[i,i],_σ[i,i],_λr[i,i],_λa[i,i],_Mwij[i,i],_σeff[i,i])
        _d = @f(d,i,_data)
    end
    _d
end

function σeff(model::SAFTVRQMieModel, V, T, z)
    _σ = model.params.sigma.values
    _λr = model.params.lambda_r.values
    _λa = model.params.lambda_a.values
    Mwij = model.params.Mw.values
    #Dij = ħ^2/(12*k_B*T*Mwij/N_A*σ[i,j]^2)

    #f(x) = Dij^2*@f(Q2,λr[i,j])+Dij*@f(Q1,λr[i,j])*x^2+x^4-x^(λr[i,j]-λa[i,j])*(Dij^2*@f(Q2,λa[i,j])+Dij*@f(Q1,λa[i,j])*x^2+x^4)
    #g(x) = 2*Dij*@f(Q1,λr[i,j])*x+4x^3-x^(λr[i,j]-λa[i,j]-1)*((λr[i,j]-λa[i,j])*Dij^2*@f(Q2,λa[i,j])+(λr[i,j]-λa[i,j]+2)*Dij*@f(Q1,λa[i,j])*x^2+(λr[i,j]-λa[i,j]+4)*x^4)
    #h(x) = 2*Dij*@f(Q1,λr[i,j])+12x^2-x^(λr[i,j]-λa[i,j]-2)*((λr[i,j]-λa[i,j]-1)*(λr[i,j]-λa[i,j])*Dij^2*@f(Q2,λa[i,j])+((λr[i,j]-λa[i,j]-1)*(λr[i,j]-λa[i,j]+2)+2*(λr[i,j]-λa[i,j]+2))*Dij*@f(Q1,λa[i,j])*x^2+((λr[i,j]-λa[i,j]-1)*(λr[i,j]-λa[i,j]+4)+4*(λr[i,j]-λa[i,j]+4))*x^4)

    function fgh(x,λa,λr,σ,T,Mw)
        Q1λr = @f(Q1,λr)
        Q1λa = @f(Q1,λa)
        Q2λr = @f(Q2,λr)
        Q2λa = @f(Q2,λa)
        Dij = ħ^2/(12*k_B*T*Mw/N_A*σ^2)
        xΔλ =x^(λr-λa)
        x2 = x*x
        x4 = x2*x2
        f = Dij^2*Q2λr+Dij*Q1λr*x2+x4-xΔλ*(Dij^2*Q2λa+Dij*Q1λa*x2+x4)
        g = 2*Dij*Q1λr*x+4x^3-(xΔλ/x)*((λr-λa)*Dij^2*Q2λa+(λr-λa+2)*Dij*Q1λa*x2+(λr-λa+4)*x4)
        h = 2*Dij*Q1λr+12x^2-x^(λr-λa-2)*((λr-λa-1)*(λr-λa)*Dij^2*Q2λa+((λr-λa-1)*(λr-λa+2)+2*(λr-λa+2))*Dij*Q1λa*x2+((λr-λa-1)*(λr-λa+4)+4*(λr-λa+4))*x4)
        return f,g,h
    end
    
    _σeff = zeros(typeof(T),size(_σ))
    n1,n2 = size(_σ)
    for i in 1:n1
        f0 = x -> fgh(x,_λa[i,i],_λr[i,i],_σ[i,i],T,Mwij[i,i])
        _σeff[i,i] =  _σ[i,i]*Solvers.halley(f0,one(T))
        for j in 1:i-1
        f0 = x -> fgh(x,_λa[i,j],_λr[i,j],_σ[i,j],T,Mwij[i,j])
        _σeff[i,j] =  _σ[i,j]*Solvers.halley(f0,one(T))
        _σeff[j,i] = _σeff[i,j]
        end
    end
    return _σeff
end

function ϵeff(model::SAFTVRQMieModel, V, T, z)
    ϵ = model.params.epsilon.values
    _σ = model.params.sigma.values
    _λr = model.params.lambda_r.values
    _λa = model.params.lambda_a.values
    Mwij = model.params.Mw.values#[i] + model.params.Mw.values[j])/2 # check
    #Dij = ħ^2/(12*k_B*T*Mwij/N_A*σ[i,j]^2)

    function fgh(x,λa,λr,σ,T,Mw)
        Q1λr = @f(Q1,λr)
        Q1λa = @f(Q1,λa)
        Q2λr = @f(Q2,λr)
        Q2λa = @f(Q2,λa)
        Dij = ħ^2/(12*k_B*T*Mw/N_A*σ^2)
        
        xΔλ =x^(λr-λa)
        x2 = x*x
        x4 = x2*x2
        f = (λr+4)*Dij^2*Q2λr+(λr+2)*Dij*Q1λr*x2+λr*x4-x^(λr-λa)*((λa+4)*Dij^2*Q2λa+(λa+2)*Dij*Q1λa*x2+λa*x4)
        g = 2*(λr+2)*Dij*Q1λr*x+4λr*x^3-x^(λr-λa-1)*((λa+4)*(λr-λa)*Dij^2*Q2λa+(λa+2)*(λr-λa+2)*Dij*Q1λa*x2+λa*(λr-λa+4)*x4)
        h = 2*(λr+2)*Dij*Q1λr+12λr*x2-x^(λr-λa-2)*((λa+4)*(λr-λa-1)*(λr-λa)*Dij^2*Q2λa+(λa+2)*((λr-λa-1)*(λr-λa+2)+2*(λr-λa+2))*Dij*Q1λa*x2+λa*((λr-λa-1)*(λr-λa+4)+4*(λr-λa+4))*x4)
    return f,g,h
    end

    function u(x,λa,λr,σ,T,Mw)
        Q1λr = @f(Q1,λr)
        Q1λa = @f(Q1,λa)
        Q2λr = @f(Q2,λr)
        Q2λa = @f(Q2,λa)
        _C = @f(Cλ,λa,λr)
        Dij = ħ^2/(12*k_B*T*Mw/N_A*σ^2)
        return _C*((x^-λr-x^-λa)+
        Dij*(Q1λr*x^-(λr+2)-Q1λa*x^-(λa+2))+
        Dij^2*(Q2λr*x^-(λr+4)-Q2λa*x^-(λa+4)))
    end
    #f(x) = (λr[i,j]+4)*Dij^2*@f(Q2,λr[i,j])+(λr[i,j]+2)*Dij*@f(Q1,λr[i,j])*x^2+λr[i,j]*x^4-x^(λr[i,j]-λa[i,j])*((λa[i,j]+4)*Dij^2*@f(Q2,λa[i,j])+(λa[i,j]+2)*Dij*@f(Q1,λa[i,j])*x^2+λa[i,j]*x^4)
    #g(x) = 2*(λr[i,j]+2)*Dij*@f(Q1,λr[i,j])*x+4λr[i,j]*x^3-x^(λr[i,j]-λa[i,j]-1)*((λa[i,j]+4)*(λr[i,j]-λa[i,j])*Dij^2*@f(Q2,λa[i,j])+(λa[i,j]+2)*(λr[i,j]-λa[i,j]+2)*Dij*@f(Q1,λa[i,j])*x^2+λa[i,j]*(λr[i,j]-λa[i,j]+4)*x^4)
    #h(x) = 2*(λr[i,j]+2)*Dij*@f(Q1,λr[i,j])+12λr[i,j]*x^2-x^(λr[i,j]-λa[i,j]-2)*((λa[i,j]+4)*(λr[i,j]-λa[i,j]-1)*(λr[i,j]-λa[i,j])*Dij^2*@f(Q2,λa[i,j])+(λa[i,j]+2)*((λr[i,j]-λa[i,j]-1)*(λr[i,j]-λa[i,j]+2)+2*(λr[i,j]-λa[i,j]+2))*Dij*@f(Q1,λa[i,j])*x^2+λa[i,j]*((λr[i,j]-λa[i,j]-1)*(λr[i,j]-λa[i,j]+4)+4*(λr[i,j]-λa[i,j]+4))*x^4)

   # σ_min = @f(Halley,f,g,h,(λr[i,j]/λa[i,j])^(1/(λr[i,j]-λa[i,j])))
    #return -ϵ[i,j]*u(σ_min)

    _ϵeff = zeros(typeof(T),size(_σ))
    n1,n2 = size(_σ)
    for i in 1:n1
        f0 = x -> fgh(x,_λa[i,i],_λr[i,i],_σ[i,i],T,Mwij[i,i])
        x0 = (_λr[i,i]/_λa[i,i])^(1/(_λr[i,i]-_λa[i,i]))
        _σmin = Solvers.halley(f0,one(T)*x0)
        uij =u(_σmin,_λa[i,i],_λr[i,i],_σ[i,i],T,Mwij[i,i])
        _ϵeff[i,i] =  -ϵ[i,i]*uij
        for j in 1:i-1
            f0 = x -> fgh(x,_λa[i,j],_λr[i,j],_σ[i,j],T,Mwij[i,j])
            x0 = (_λr[i,j]/_λa[i,j])^(1/(_λr[i,j]-_λa[i,j]))
            _σmin = Solvers.halley(f0,one(T)*x0)
            uij =u(_σmin,_λa[i,j],_λr[i,j],_σ[i,j],T,_Mw[i,j])
            _ϵeff[i,j] =  -ϵ[i,j]*uij
            _ϵeff[j,i] = _ϵeff[i,j]
        end
    end
    return _ϵeff
end

#=
function a_1(model::SAFTVRQMieModel, V, T, z, i, j,_data = @f(data))
    _σeff,_ϵff,_d,_ρ_S,ζi,_ζ_X,_ζst = _data
    ϵ = model.params.epsilon.values
    σ = model.params.sigma.values
    λr = model.params.lambda_r.values
    λa = model.params.lambda_a.values
    Mwij = (model.params.Mw.values[i] + model.params.Mw.values[j])/2 # check
    Dij = ħ^2/(12*k_B*T*Mwij/N_A*σ[i,j]^2)
    x_0ij = @f(x_0,i,j)
    dij = 0.5*(_d[i] + _d[j])
    x_0effij = _σeff[i,j]/dij
    return 2*π*ϵ[i,j]*dij^3*@f(Cλ,λa[i,j],λr[i,j])*@f(ρ_S)*
    ( (x_0ij^λa[i,j]*(@f(aS_1,λa[i,j])+@f(B,λa[i,j],x_0effij))-
        x_0ij^λr[i,j]*(@f(aS_1,λr[i,j])+@f(B,λr[i,j],x_0effij)))+
        (x_0ij^(λa[i,j]+2)*@f(Q1,λa[i,j])*(@f(aS_1,λa[i,j]+2)+@f(B,λa[i,j]+2,x_0effij))-
         x_0ij^(λr[i,j]+2)*@f(Q1,λr[i,j])*(@f(aS_1,λr[i,j]+2)+@f(B,λr[i,j]+2,x_0effij)))*Dij+
        (x_0ij^(λa[i,j]+4)*@f(Q2,λa[i,j])*(@f(aS_1,λa[i,j]+4)+@f(B,λa[i,j]+4,x_0effij))-
         x_0ij^(λr[i,j]+4)*@f(Q2,λr[i,j])*(@f(aS_1,λr[i,j]+4)+@f(B,λr[i,j]+4,x_0effij)))*Dij^2 )
end
=#

#=
function x_0eff(model::SAFTVRQMieModel, V, T, z, i, j)
    return @f(σeff,i,j)/@f(d,i,j)
end

function a_2(model::SAFTVRQMieModel, V, T, z, i, j,_data = @f(data))
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



function a_3(model::SAFTVRQMieModel, V, T, z, i, j)
    λr = model.params.lambda_r.values
    λa = model.params.lambda_a.values
    σ = model.params.sigma.values
    ϵ = model.params.epsilon.values
    σeffij = @f(σeff,i,j)
    ϵeffij = @f(ϵeff,i,j)
    Mwij = model.params.Mw.values[i,j]
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
=#
function a_disp(model::SAFTVRQMieModel, V, T, z,_data = @f(data))
    _σeff,_ϵff,vrdata= _data
    _d,_ρ_S,ζi,_ζ_X,_ζst  = vrdata
    comps = @comps
    l = length(comps)
    ∑z = ∑(z)
    m = model.params.segment.values
    _ϵ = model.params.epsilon.values
    _λr = model.params.lambda_r.values
    _λa = model.params.lambda_a.values
    _σ = model.params.sigma.values
    Mw = model.params.Mw.values
    m̄ = dot(z, m)
    m̄inv = 1/m̄
    a₁ = zero(V+T+first(z))
    a₂ = a₁
    a₃ = a₁
    _ζst5 = _ζst^5
    _ζst8 = _ζst^8
    _KHS = @f(KHS,_ζ_X,_ρ_S)
    for i ∈ comps
        j = i
        x_Si = z[i]*m[i]*m̄inv
        x_Sj = x_Si
        ϵ,λa,λr,σ,di= _ϵ[i,j],_λa[i,i],_λr[i,i],_σ[i,i],_d[i]
        σeff,ϵff = _σeff[i,j],_ϵff[i,j]
        _C = @f(Cλ,λa,λr)
        dij3 = di^3
        x_0ij = σ/di
        x_0effij = σeff/di
        Mwij = Mw[i,j]
        Dij = ħ^2/(12*k_B*T*Mwij/N_A*σ^2)
        Q1λr = @f(Q1,λr)
        Q1λa = @f(Q1,λa)
        Q2λr = @f(Q2,λr)
        Q2λa = @f(Q2,λa)
        
        #calculations for a1 - diagonal
        a1_ij = 2*π*ϵ*dij3*_C*_ρ_S*
            ( (x_0ij^λa*(@f(aS_1,λa,_ζ_X)+@f(B,λa,x_0effij,_ζ_X))-
            x_0ij^λr*(@f(aS_1,λr,_ζ_X)+@f(B,λr,x_0effij,_ζ_X)))+
            (x_0ij^(λa+2)*Q1λa*(@f(aS_1,λa+2,_ζ_X)+@f(B,λa+2,x_0effij,_ζ_X))-
            x_0ij^(λr+2)*Q1λr*(@f(aS_1,λr+2,_ζ_X)+@f(B,λr+2,x_0effij,_ζ_X)))*Dij+
            (x_0ij^(λa+4)*Q2λa*(@f(aS_1,λa+4,_ζ_X)+@f(B,λa+4,x_0effij,_ζ_X))-
            x_0ij^(λr+4)*Q2λr*(@f(aS_1,λr+4,_ζ_X)+@f(B,λr+4,x_0effij,_ζ_X)))*Dij^2 )
        #calculations for a2 - diagonal
        σcoeff =σ/σeff 
        α = _C*ϵ/ϵff*
            ((σcoeff^λa/(λa-3)-σcoeff^λr/(λr-3))+
            Dij*(σcoeff^(2+λa)*Q1λa/(λa-1)-
            σcoeff^(2+λr)*Q1λr/(λr-1))+
            Dij^2*(σcoeff^(4+λa)*Q2λa/(λa+1)-
            σcoeff^(4+λr)*Q2λr/(λr+1)))
        f1,f2,f3,f4,f5,f6 = @f(f123456,α)
         _χ = f1*_ζst+f2*_ζst5+f3*_ζst8
        a2_ij = π*_KHS*(1+_χ)*_ρ_S*ϵ^2*dij3*_C^2*(x_0ij^(2*λa)*(@f(aS_1,2*λa)+@f(B,2*λa,x_0effij,_ζ_X))-
            x_0ij^(λa+λr)*2*(@f(aS_1,λa+λr)+@f(B,λa+λr,x_0effij,_ζ_X))+
            x_0ij^(2*λr)*(@f(aS_1,2*λr)+@f(B,2*λr,x_0effij,_ζ_X))+
            x_0ij^(2*λa+2)*2*Q1λa*(@f(aS_1,2*λa+2)+@f(B,2*λa+2,x_0effij,_ζ_X))*Dij+
            x_0ij^(2*λr+2)*2*Q1λr*(@f(aS_1,2*λr+2)+@f(B,2*λr+2,x_0effij,_ζ_X))*Dij-
            x_0ij^(λa+λr+2)*2*(Q1λa+Q1λr)*(@f(aS_1,λa+λr+2)+@f(B,λa+λr+2,x_0effij,_ζ_X))*Dij+
            x_0ij^(2*λa+4)*Q1λa^2*(@f(aS_1,2*λa+4,_ζ_X)+@f(B,2*λa+4,x_0effij,_ζ_X))*Dij^2+
            x_0ij^(2*λr+4)*Q1λr^2*(@f(aS_1,2*λr+4,_ζ_X)+@f(B,2*λr+4,x_0effij,_ζ_X))*Dij^2-
            x_0ij^(λa+λr+4)*(2*Q1λa*Q1λr)*(@f(aS_1,λa+λr+4,_ζ_X)+@f(B,λa+λr+4,x_0effij,_ζ_X))*Dij^2+
            x_0ij^(2*λa+4)*2*Q2λa*(@f(aS_1,2*λa+4,_ζ_X)+@f(B,2*λa+4,x_0effij,_ζ_X))*Dij^2+
            x_0ij^(2*λr+4)*2*Q2λr*(@f(aS_1,2*λr+4,_ζ_X)+@f(B,2*λr+4,x_0effij,_ζ_X))*Dij^2-
            x_0ij^(λa+λr+4)*2*(Q2λa+Q2λr)*(@f(aS_1,λa+λr+4,_ζ_X)+@f(B,λa+λr+4,x_0effij,_ζ_X))*Dij^2+
            x_0ij^(2*λa+6)*2*(Q1λa*Q2λa)*(@f(aS_1,2*λa+6)+@f(B,2*λa+6,x_0effij,_ζ_X))*Dij^3+
            x_0ij^(2*λr+6)*2*(Q1λr*Q2λr)*(@f(aS_1,2*λr+6)+@f(B,2*λr+6,x_0effij,_ζ_X))*Dij^3-
            x_0ij^(λa+λr+6)*2*(Q1λr*Q2λa+Q1λa*Q2λr)*(@f(aS_1,λa+λr+6)+@f(B,λa+λr+6,x_0effij,_ζ_X))*Dij^3+
            x_0ij^(2*λa+8)*Q2λa^2*(@f(aS_1,2*λa+8,_ζ_X)+@f(B,2*λa+8,x_0effij,_ζ_X))*Dij^4+
            x_0ij^(2*λr+8)*Q2λr^2*(@f(aS_1,2*λr+8,_ζ_X)+@f(B,2*λr+8,x_0effij,_ζ_X))*Dij^4-
            x_0ij^(λa+λr+8)*(2*Q2λa*Q2λr)*(@f(aS_1,λa+λr+8,_ζ_X)+@f(B,λa+λr+8,x_0effij,_ζ_X))*Dij^4)
    
        
        #calculations for a3 - diagonal
        a3_ij = -ϵff^3*f4*_ζst * exp(f5*_ζst+f6*_ζst^2)
        #adding - diagonal
        a₁ += a1_ij*x_Si*x_Si
        a₂ += a2_ij*x_Si*x_Si
        a₃ += a3_ij*x_Si*x_Si
        for j ∈ (i+1):l
            x_Sj = z[j]*m[j]*m̄inv   
            ϵ,λa,λr,σ,dj= _ϵ[i,j],_λa[i,j],_λr[i,j],_σ[i,j],_d[j]
            σeff,ϵff = _σeff[i,j],_ϵff[i,j]
            _C = @f(Cλ,λa,λr)
            dij = 0.5*(di+dj)
            dij3 = dij^3
            x_0ij = σ/dij
            x_0effij = σeff/dij
            Q1λr = @f(Q1,λr)
            Q1λa = @f(Q1,λa)
            Q2λr = @f(Q2,λr)
            Q2λa = @f(Q2,λa)
            Mwij = Mw[i,j]
            Dij = ħ^2/(12*k_B*T*Mwij/N_A*σ^2)
    
            #calculations for a1
            a1_ij = 2*π*ϵ*dij3*_C*ρ_S*
            ( (x_0ij^λa*(@f(aS_1,λa,_ζ_X)+@f(B,λa,x_0effij,_ζ_X))-
                x_0ij^λr*(@f(aS_1,λr,_ζ_X)+@f(B,λr,x_0effij,_ζ_X)))+
                (x_0ij^(λa+2)*Q1λa*(@f(aS_1,λa+2,_ζ_X)+@f(B,λa+2,x_0effij,_ζ_X))-
                 x_0ij^(λr+2)*Q1λr*(@f(aS_1,λr+2,_ζ_X)+@f(B,λr+2,x_0effij,_ζ_X)))*Dij+
                (x_0ij^(λa+4)*Q2λa*(@f(aS_1,λa+4,_ζ_X)+@f(B,λa+4,x_0effij,_ζ_X))-
                 x_0ij^(λr+4)*Q2λr*(@f(aS_1,λr+4,_ζ_X)+@f(B,λr+4,x_0effij,_ζ_X)))*Dij^2 )
    
            #calculations for a2
            σcoeff =σ/σeff 
            α = _C*ϵ/ϵff*
                ((σcoeff^λa/(λa-3)-σcoeff^λr/(λr-3))+
                Dij*(σcoeff^(2+λa)*Q1λa/(λa-1)-
                σcoeff^(2+λr)*Q1λr/(λr-1))+
                Dij^2*(σcoeff^(4+λa)*Q2λa/(λa+1)-
                σcoeff^(4+λr)*Q2λr/(λr+1)))
            f1,f2,f3,f4,f5,f6 = @f(f123456,α)
             _χ = f1*_ζst+f2*_ζst5+f3*_ζst8
             a2_ij = π*_KHS*(1+_χ)*_ρ_S*ϵ^2*dij3*_C^2*(x_0ij^(2*λa)*(@f(aS_1,2*λa)+@f(B,2*λa,x_0effij,_ζ_X))-
             x_0ij^(λa+λr)*2*(@f(aS_1,λa+λr)+@f(B,λa+λr,x_0effij,_ζ_X))+
             x_0ij^(2*λr)*(@f(aS_1,2*λr)+@f(B,2*λr,x_0effij,_ζ_X))+
             x_0ij^(2*λa+2)*2*Q1λa*(@f(aS_1,2*λa+2)+@f(B,2*λa+2,x_0effij,_ζ_X))*Dij+
             x_0ij^(2*λr+2)*2*Q1λr*(@f(aS_1,2*λr+2)+@f(B,2*λr+2,x_0effij,_ζ_X))*Dij-
             x_0ij^(λa+λr+2)*2*(Q1λa+Q1λr)*(@f(aS_1,λa+λr+2)+@f(B,λa+λr+2,x_0effij,_ζ_X))*Dij+
             x_0ij^(2*λa+4)*Q1λa^2*(@f(aS_1,2*λa+4,_ζ_X)+@f(B,2*λa+4,x_0effij,_ζ_X))*Dij^2+
             x_0ij^(2*λr+4)*Q1λr^2*(@f(aS_1,2*λr+4,_ζ_X)+@f(B,2*λr+4,x_0effij,_ζ_X))*Dij^2-
             x_0ij^(λa+λr+4)*(2*Q1λa*Q1λr)*(@f(aS_1,λa+λr+4,_ζ_X)+@f(B,λa+λr+4,x_0effij,_ζ_X))*Dij^2+
             x_0ij^(2*λa+4)*2*Q2λa*(@f(aS_1,2*λa+4,_ζ_X)+@f(B,2*λa+4,x_0effij,_ζ_X))*Dij^2+
             x_0ij^(2*λr+4)*2*Q2λr*(@f(aS_1,2*λr+4,_ζ_X)+@f(B,2*λr+4,x_0effij,_ζ_X))*Dij^2-
             x_0ij^(λa+λr+4)*2*(Q2λa+Q2λr)*(@f(aS_1,λa+λr+4,_ζ_X)+@f(B,λa+λr+4,x_0effij,_ζ_X))*Dij^2+
             x_0ij^(2*λa+6)*2*(Q1λa*Q2λa)*(@f(aS_1,2*λa+6)+@f(B,2*λa+6,x_0effij,_ζ_X))*Dij^3+
             x_0ij^(2*λr+6)*2*(Q1λr*Q2λr)*(@f(aS_1,2*λr+6)+@f(B,2*λr+6,x_0effij,_ζ_X))*Dij^3-
             x_0ij^(λa+λr+6)*2*(Q1λr*Q2λa+Q1λa*Q2λr)*(@f(aS_1,λa+λr+6)+@f(B,λa+λr+6,x_0effij,_ζ_X))*Dij^3+
             x_0ij^(2*λa+8)*Q2λa^2*(@f(aS_1,2*λa+8,_ζ_X)+@f(B,2*λa+8,x_0effij,_ζ_X))*Dij^4+
             x_0ij^(2*λr+8)*Q2λr^2*(@f(aS_1,2*λr+8,_ζ_X)+@f(B,2*λr+8,x_0effij,_ζ_X))*Dij^4-
             x_0ij^(λa+λr+8)*(2*Q2λa*Q2λr)*(@f(aS_1,λa+λr+8,_ζ_X)+@f(B,λa+λr+8,x_0effij,_ζ_X))*Dij^4)
    
            #calculations for a3
            a3_ij = -ϵff^3*f4*_ζst * exp(f5*_ζst+f6*_ζst^2)
            #adding
            a₁ += 2*a1_ij*x_Si*x_Sj
            a₂ += 2*a2_ij*x_Si*x_Sj
            a₃ += 2*a3_ij*x_Si*x_Sj            
        end
    end
    a₁ = a₁*m̄/T/∑z
    a₂ = a₂*m̄/(T*T)/∑z
    a₃ = a₃*m̄/(T*T*T)/∑z
    adisp =  a₁ + a₂ + a₃ 
    return adisp
end

const SAFTVRQMieconsts = (
    x = [-0.9956571630258080807355,-0.973906528517171720078,-0.9301574913557082260012,-0.8650633666889845107321,-0.7808177265864168970637,-0.6794095682990244062343,-0.562757134668604683339,
        -0.4333953941292471907993,-0.294392862701460198131,-0.1488743389816312108848,0,0.1488743389816312108848,0.2943928627014601981311,0.4333953941292471907993,
        0.562757134668604683339,0.6794095682990244062343,0.7808177265864168970637,0.865063366688984510732,0.9301574913557082260012,0.973906528517171720078,0.9956571630258080807355],
    w = [0.0116946388673718742781,0.0325581623079647274788,0.0547558965743519960314,0.075039674810919952767,0.093125454583697605535,0.1093871588022976418992,0.123491976262065851078,
        0.134709217311473325928,0.142775938577060080797,0.1477391049013384913748,0.149445554002916905665,0.1477391049013384913748,0.1427759385770600807971,0.134709217311473325928,
        0.123491976262065851078,0.109387158802297641899,0.093125454583697605535,0.075039674810919952767,0.05475589657435199603138,0.032558162307964727479,0.0116946388673718742781],
)
