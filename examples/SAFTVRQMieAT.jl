struct SAFTVRQMieATParam <: EoSParam
    segment::SingleParam{Float64}
    sigma::PairParam{Float64}
    lambda_a::PairParam{Float64}
    lambda_r::PairParam{Float64}
    epsilon::PairParam{Float64}
    Mw::SingleParam{Float64}
    nu::SingleParam{Float64}
end

abstract type SAFTVRQMieATModel <: SAFTVRMieModel end

@newmodel SAFTVRQMieAT SAFTVRQMieATModel SAFTVRQMieATParam
export SAFTVRQMieAT

# Although optional, there should be some additional argument if the user would like to modify certain aspects of the model, such as the ideal model used
function SAFTVRQMieAT(components; 
    idealmodel=BasicIdeal, 
    userlocations=String[], 
    ideal_userlocations=String[], 
    verbose=false)

    # The next step is to collect all the parameters from the csv files. The getparams function can handle this. 
    # It is also useful to allow for custom userlocations to be fed into this function
    params = Clapeyron.getparams(components, ["properties/molarmass.csv"]; userlocations=append!(userlocations,ideal_userlocations), verbose=verbose)

    # The next step is to transform the parameters from the database into the form in which they are useable within the equation of state (e.g. applying combining rules)
    params["Mw"].values .*= 1E-3
    Mw = params["Mw"]
    segment = params["m"]
    nu = params["nu"]
    params["sigma"].values .*= 1E-10
    sigma = Clapeyron.sigma_LorentzBerthelot(params["sigma"])
    epsilon = Clapeyron.epsilon_HudsenMcCoubrey(params["epsilon"], sigma)
    lambda_a = Clapeyron.lambda_LorentzBerthelot(params["lambda_a"])
    lambda_r = Clapeyron.lambda_LorentzBerthelot(params["lambda_r"])

    # The struct with all of our parameters is then formed
    packagedparams = SAFTVRQMieATParam(segment, sigma, lambda_a, lambda_r, epsilon, Mw, nu)

    # If the model is based on an existing one, it's generally a good idea to add the reference DOI here
    references = [""]

    # We then construct and return the model
    model = SAFTVRQMieAT(packagedparams, idealmodel;ideal_userlocations=ideal_userlocations, references=references, verbose=verbose)
    return model
end

function d(model::SAFTVRQMieATModel, V, T, z)
    ρ  = N_A*sum(z)/V
    λ  = 0.85
    ν  = model.params.nu.values
    ϵ = model.params.epsilon.diagvalues
    σ = model.params.sigma.diagvalues
    λr = model.params.lambda_r.diagvalues
    λa = model.params.lambda_a.diagvalues
    return d.(model,V,T,Ref(z),λa,λr,ϵ,σ,ν,ρ,λ)
end

function d(model::SAFTVRQMieATModel, V, T, z, λa,λr,ϵ,σ,ν,ρ,λ)
    u = Clapeyron.SAFTVRMieconsts.u
    w = Clapeyron.SAFTVRMieconsts.w
    θ = @f(Clapeyron.Cλ,λa,λr)*ϵ/T*(1-ρ*λ*ν*σ^3)
    σ*(1-∑(w[j]*(θ/(θ+u[j]))^(1/λr)*(exp(θ*(1/(θ/(θ+u[j]))^(λa/λr)-1))/(u[j]+θ)/λr) for j ∈ 1:5))
end

function a_disp(model::SAFTVRQMieATModel, V, T, z,_data = @f(data))
    _d,ρS,ζi,_ζ_X,_ζst,_ = _data
    comps = @Clapeyron.comps
    l = length(comps)
    ∑z = ∑(z)
    m = model.params.segment.values
    _ϵ = model.params.epsilon.values
    _λr = model.params.lambda_r.values
    _λa = model.params.lambda_a.values
    _σ = model.params.sigma.values
    _ν = model.params.nu.values
    m̄ = Clapeyron.dot(z, m)
    m̄inv = 1/m̄
    a₁ = zero(V+T+first(z))
    a₂ = a₁
    a₃ = a₁
    _ζst5 = _ζst^5
    _ζst8 = _ζst^8
    _KHS = @f(KHS,_ζ_X,ρS)
    
    ρ  = N_A*sum(z)/V
    λ  = 0.85
    for i ∈ comps
        j = i
        x_Si = z[i]*m[i]*m̄inv
        x_Sj = x_Si
        ϵ = _ϵ[i,j]
        λa = _λa[i,i]
        λr = _λr[i,i] 
        σ = _σ[i,i]
        ν  = _ν[i]
        AT_corr = (1-ρ*λ*ν*σ^3)
        _C = @f(Cλ,λa,λr)
        dij = _d[i]
        dij3 = dij^3
        x_0ij = σ/dij
        #calculations for a1 - diagonal
        aS_1_a = @f(aS_1,λa,_ζ_X)
        aS_1_r = @f(aS_1,λr,_ζ_X)
        B_a = @f(B,λa,x_0ij,_ζ_X)
        B_r = @f(B,λr,x_0ij,_ζ_X)
        a1_ij = (2*π*AT_corr*ϵ*dij3)*_C*ρS*
        (x_0ij^λa*(aS_1_a+B_a) - x_0ij^λr*(aS_1_r+B_r))

        #calculations for a2 - diagonal
        aS_1_2a = @f(aS_1,2*λa,_ζ_X)
        aS_1_2r = @f(aS_1,2*λr,_ζ_X)
        aS_1_ar = @f(aS_1,λa+λr,_ζ_X)
        B_2a = @f(B,2*λa,x_0ij,_ζ_X)
        B_2r = @f(B,2*λr,x_0ij,_ζ_X)
        B_ar = @f(B,λr+λa,x_0ij,_ζ_X)
        α = _C*(1/(λa-3)-1/(λr-3))
        f1,f2,f3,f4,f5,f6 = @f(f123456,α)
        _χ = f1*_ζst+f2*_ζst5+f3*_ζst8
        a2_ij = π*_KHS*(1+_χ)*ρS*AT_corr^2*ϵ^2*dij3*_C^2 *
        (x_0ij^(2*λa)*(aS_1_2a+B_2a)
        - 2*x_0ij^(λa+λr)*(aS_1_ar+B_ar)
        + x_0ij^(2*λr)*(aS_1_2r+B_2r))
        
        #calculations for a3 - diagonal
        a3_ij = -ϵ^3*AT_corr^3*f4*_ζst * exp(f5*_ζst+f6*_ζst^2)
        #adding - diagonal
        a₁ += a1_ij*x_Si*x_Si
        a₂ += a2_ij*x_Si*x_Si
        a₃ += a3_ij*x_Si*x_Si
    end
    a₁ = a₁*m̄/T/∑z
    a₂ = a₂*m̄/(T*T)/∑z
    a₃ = a₃*m̄/(T*T*T)/∑z
    #@show (a₁,a₂,a₃)
    adisp =  a₁ + a₂ + a₃ 
    return adisp
end

function a_res(model::SAFTVRQMieATModel, V, T, z)
    _data = @f(data)
    return @f(a_hs,_data)+@f(a_disp,_data)
end