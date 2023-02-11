struct SAFTVRMieParam <: EoSParam
    Mw::SingleParam{Float64}
    segment::SingleParam{Float64}
    sigma::PairParam{Float64}
    lambda_a::PairParam{Float64}
    lambda_r::PairParam{Float64}
    epsilon::PairParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
end

abstract type SAFTVRMieModel <: SAFTModel end
@newmodel SAFTVRMie SAFTVRMieModel SAFTVRMieParam

"""
    SAFTVRMieModel <: SAFTModel

    SAFTVRMie(components; 
    idealmodel=BasicIdeal,
    userlocations=String[],
    ideal_userlocations=String[],
    verbose=false,
    assoc_options = AssocOptions())

## Input parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Single Parameter (`Float64`) - Segment Diameter [`A°`]
- `epsilon`: Single Parameter (`Float64`) - Reduced dispersion energy  `[K]`
- `lambda_a`: Pair Parameter (`Float64`) - Atractive range parameter (no units)
- `lambda_r`: Pair Parameter (`Float64`) - Repulsive range parameter (no units)
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Paramater (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume

## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Pair Parameter (`Float64`) - Mixed segment Diameter `[m]`
- `lambda_a`: Pair Parameter (`Float64`) - Atractive range parameter (no units)
- `lambda_r`: Pair Parameter (`Float64`) - Repulsive range parameter (no units)
- `epsilon`: Pair Parameter (`Float64`) - Mixed reduced dispersion energy`[K]`
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume

## Input models
- `idealmodel`: Ideal Model

## Description

SAFT-VR with Mie potential

## References
1. Lafitte, T., Apostolakou, A., Avendaño, C., Galindo, A., Adjiman, C. S., Müller, E. A., & Jackson, G. (2013). Accurate statistical associating fluid theory for chain molecules formed from Mie segments. The Journal of Chemical Physics, 139(15), 154504. [doi:10.1063/1.4819786](https://doi.org/10.1063/1.4819786)
2. Dufal, S., Lafitte, T., Haslam, A. J., Galindo, A., Clark, G. N. I., Vega, C., & Jackson, G. (2015). The A in SAFT: developing the contribution of association to the Helmholtz free energy within a Wertheim TPT1 treatment of generic Mie fluids. Molecular Physics, 113(9–10), 948–984. [doi:10.1080/00268976.2015.1029027](https://doi.org/10.1080/00268976.2015.1029027)
"""
SAFTVRMie

export SAFTVRMie

function SAFTVRMie(components;
    idealmodel=BasicIdeal,
    userlocations=String[],
    ideal_userlocations=String[],
    verbose=false,
    assoc_options = AssocOptions(), kwargs...)
    params,sites = getparams(components, ["SAFT/SAFTVRMie", "properties/molarmass.csv"]; userlocations=userlocations, verbose=verbose)

    Mw = params["Mw"]
    segment = params["segment"]
    params["sigma"].values .*= 1E-10
    sigma = sigma_LorentzBerthelot(params["sigma"])
    epsilon = epsilon_HudsenMcCoubrey(params["epsilon"], sigma)
    lambda_a = lambda_LorentzBerthelot(params["lambda_a"])
    lambda_r = lambda_LorentzBerthelot(params["lambda_r"])
    epsilon_assoc = params["epsilon_assoc"]
    bondvol = params["bondvol"]
    bondvol,epsilon_assoc = assoc_mix(bondvol,epsilon_assoc,sigma,assoc_options) #combining rules for association

    packagedparams = SAFTVRMieParam(Mw, segment, sigma, lambda_a, lambda_r, epsilon, epsilon_assoc, bondvol)
    references = ["10.1063/1.4819786", "10.1080/00268976.2015.1029027"]

    model = SAFTVRMie(packagedparams, sites, idealmodel; ideal_userlocations, references, verbose, assoc_options)
    return model
end

function recombine_impl!(model::SAFTVRMieModel)
    assoc_options = model.assoc_options
    sigma = model.params.sigma
    epsilon = model.params.epsilon
    lambda_a = model.params.lambda_a
    lambda_r = model.params.lambda_r
    
    epsilon_assoc = model.params.epsilon_assoc
    bondvol = model.params.bondvol
    bondvol,epsilon_assoc = assoc_mix(bondvol,epsilon_assoc,sigma,assoc_options) #combining rules for association

    model.params.epsilon_assoc.values.values[:] = epsilon_assoc.values.values
    model.params.bondvol.values.values[:] = bondvol.values.values

    sigma = sigma_LorentzBerthelot!(sigma)
    epsilon = epsilon_HudsenMcCoubrey!(epsilon,sigma)
    lambda_a = lambda_LorentzBerthelot!(lambda_a)
    lambda_r = lambda_LorentzBerthelot!(lambda_r)
    return model
end

function x0_volume_liquid(model::SAFTVRMieModel,T,z)
    v_lb = lb_volume(model,z)
    return v_lb*1.5
end

function data(model::SAFTVRMieModel, V, T, z)
    m̄ = dot(z,model.params.segment)
    _d = @f(d)
    ζi = @f(ζ0123,_d)
    _ζ_X,σ3x = @f(ζ_X_σ3,_d,m̄)
    _ρ_S = @f(ρ_S,m̄)
    _ζst = σ3x*_ρ_S*π/6  
    return (_d,_ρ_S,ζi,_ζ_X,_ζst,σ3x,m̄)
end

function a_res(model::SAFTVRMieModel, V, T, z)
    _data = @f(data)
    return @f(a_hs,_data)+@f(a_disp,_data) + @f(a_chain,_data) + @f(a_assoc,_data)
end

#fused chain and disp calculation
function a_res(model::SAFTVRMie, V, T, z)
    _data = @f(data)
    return @f(a_hs,_data)+@f(a_dispchain,_data) + @f(a_assoc,_data)
end

function a_mono(model::SAFTVRMieModel, V, T, z,_data = @f(data))
    return @f(a_hs,_data)+@f(a_disp,_data)
end

function a_hs(model::SAFTVRMieModel, V, T, z,_data = @f(data))
    _,_,ζi,_,_,_,_ = _data
    ζ0,ζ1,ζ2,ζ3 = ζi
    N = N_A*∑(z)
    return 6*V/π/N*(3ζ1*ζ2/(1-ζ3) + ζ2^3/(ζ3*(1-ζ3)^2) + (ζ2^3/ζ3^2-ζ0)*log(1-ζ3))
end

function ρ_S(model::SAFTVRMieModel, V, T, z, m̄ = dot(z,model.params.segment.values))
    return N_A/V*m̄
end

function ζ0123(model::SAFTVRMieModel, V, T, z,_d=@f(d),m̄ = dot(z,model.params.segment.values))
    m = model.params.segment
    _0 = zero(V+T+first(z))
    ζ0,ζ1,ζ2,ζ3 = _0,_0,_0,_0
    for i ∈ 1:length(z)
        di =_d[i]
        xS = z[i]*m[i]/m̄
        ζ0 += xS
        ζ1 += xS*di
        ζ2 += xS*di*di
        ζ3 += xS*di*di*di
    end
    c = π/6*N_A*m̄/V
    ζ0,ζ1,ζ2,ζ3 = c*ζ0,c*ζ1,c*ζ2,c*ζ3
    return ζ0,ζ1,ζ2,ζ3 
end


function d(model::SAFTVRMieModel, V, T, z)
    _ϵ = diagvalues(model.params.epsilon)
    _σ = diagvalues(model.params.sigma)
    _λa = diagvalues(model.params.lambda_a)
    _λr = diagvalues(model.params.lambda_r)
    u = SAFTVRMieconsts.u
    w = SAFTVRMieconsts.w
    _0 = zero(T*1.0)
    n = length(z)
    _d = zeros(typeof(_0),n)
    for k ∈ 1:n
        ϵ = _ϵ[k]
        σ = _σ[k]
        λa = _λa[k]
        λr = _λr[k]
        θ = Cλ(model,V,T,z,λa,λr)*ϵ/T
        di = _0
        λrinv = 1/λr
        λaλr = λa/λr
        for j ∈ 1:5
            di += w[j]*(θ/(θ+u[j]))^λrinv*(exp(θ*(1/(θ/(θ+u[j]))^λaλr-1))/(u[j]+θ)/λr)
        end
        _d[k] = σ*(1-di)
    end
    return _d
end


function d(model::SAFTVRMieModel, V, T, z, λa,λr,ϵ,σ)
    u = SAFTVRMieconsts.u
    w = SAFTVRMieconsts.w
    θ = @f(Cλ,λa,λr)*ϵ/T
    σ*(1-∑(w[j]*(θ/(θ+u[j]))^(1/λr)*(exp(θ*(1/(θ/(θ+u[j]))^(λa/λr)-1))/(u[j]+θ)/λr) for j ∈ 1:5))
end


function Cλ(model::SAFTVRMieModel, V, T, z, λa, λr)
    return (λr/(λr-λa))*(λr/λa)^(λa/(λr-λa))
end

function ζ_X(model::SAFTVRMieModel, V, T, z,_d = @f(d))
    _ζ_X,σ3x = @f(ζ_X_σ3,_d)
    return _ζ_X
end

function ζ_X_σ3(model::SAFTVRMieModel, V, T, z,_d = @f(d),m̄ = dot(z,model.params.segment.values))
    m = model.params.segment.values
    m̄ = dot(z, m)
    m̄inv = 1/m̄
    σ = model.params.sigma
    ρS = N_A/V*m̄
    comps = 1:length(z)
    _ζ_X = zero(V+T+first(z))
    kρS = ρS* π/6/8 
    σ3_x = _ζ_X

    for i ∈ comps
        x_Si = z[i]*m[i]*m̄inv
        σ3_x += x_Si*x_Si*(σ[i,i]^3)
        di =_d[i]
        r1 = kρS*x_Si*x_Si*(2*di)^3         
        _ζ_X += r1
        for j ∈ 1:(i-1)
            x_Sj = z[j]*m[j]*m̄inv
            σ3_x += 2*x_Si*x_Sj*(σ[i,j]^3)
            dij = (di + _d[j])
            r1 = kρS*x_Si*x_Sj*dij^3         
            _ζ_X += 2*r1
        end
    end
   
    return _ζ_X,σ3_x
end

function aS_1(model::SAFTVRMieModel, V, T, z, λ,ζ_X_= @f(ζ_X))
    ζeff_ = @f(ζeff,λ,ζ_X_)
    return -1/(λ-3)*(1-ζeff_/2)/(1-ζeff_)^3
end

function ζeff(model::SAFTVRMieModel, V, T, z, λ,ζ_X_= @f(ζ_X))
    A = SAFTγMieconsts.A
    λ⁻¹ = one(λ)/λ
    return A * SA[one(λ); λ⁻¹; λ⁻¹*λ⁻¹; λ⁻¹*λ⁻¹*λ⁻¹] ⋅ SA[ζ_X_; ζ_X_^2; ζ_X_^3; ζ_X_^4]
end

function B(model::SAFTVRMieModel, V, T, z, λ, x_0,ζ_X_ = @f(ζ_X))
    I = (1-x_0^(3-λ))/(λ-3)
    J = (1-(λ-3)*x_0^(4-λ)+(λ-4)*x_0^(3-λ))/((λ-3)*(λ-4))
    return I*(1-ζ_X_/2)/(1-ζ_X_)^3-9*J*ζ_X_*(ζ_X_+1)/(2*(1-ζ_X_)^3)
end

function KHS(model::SAFTVRMieModel, V, T, z,ζ_X_ = @f(ζ_X),ρS=@f(ρ_S))
    return (1-ζ_X_)^4/(1+4ζ_X_+4ζ_X_^2-4ζ_X_^3+ζ_X_^4)
end

function f123456(model::SAFTVRMieModel, V, T, z, α)
    ϕ = SAFTVRMieconsts.ϕ
    _0 = zero(α)
    fa = (_0,_0,_0,_0,_0,_0)
    fb = (_0,_0,_0,_0,_0,_0)
    @inbounds for i ∈ 1:4
        ϕi = ϕ[i]::NTuple{6,Float64}
        ii = i-1
        αi = α^ii
        fa = fa .+ ϕi .*αi
    end
    @inbounds for i ∈ 5:7
        ϕi = ϕ[i]::NTuple{6,Float64}
        ii = i-4
        αi = α^ii
        fb = fb .+ ϕi .*αi
    end
    return  fa ./ (one(_0) .+ fb)
    #return sum(ϕ[i+1][m]*α^i for i ∈ 0:3)/(1+∑(ϕ[i+1][m]*α^(i-3) for i ∈ 4:6))
end

function ζst(model::SAFTVRMieModel, V, T, z,_σ = model.params.sigma)
    m = model.params.segment.values
    m̄ = dot(z, m)
    m̄inv = 1/m̄
    ρS = N_A/V*m̄
    comps = @comps
    _ζst = zero(V+T+first(z))
    for i ∈ comps
        x_Si = z[i]*m[i]*m̄inv
        _ζst += x_Si*x_Si*(_σ[i,i]^3)
        for j ∈ 1:i-1
            x_Sj = z[j]*m[j]*m̄inv
            _ζst += 2*x_Si*x_Sj*(_σ[i,j]^3)           
        end
    end

    #return π/6*@f(ρ_S)*∑(@f(x_S,i)*@f(x_S,j)*(@f(d,i)+@f(d,j))^3/8 for i ∈ comps for j ∈ comps)
    return _ζst*ρS* π/6
end

function g_HS(model::SAFTVRMieModel, V, T, z, x_0ij,ζ_X_ = @f(ζ_X))
    ζX3 = (1-ζ_X_)^3
    #evalpoly(ζ_X_,(0,42,-39,9,-2)) = (42ζ_X_-39ζ_X_^2+9ζ_X_^3-2ζ_X_^4)
    k_0 = -log(1-ζ_X_)+evalpoly(ζ_X_,(0,42,-39,9,-2))/(6*ζX3)
    #evalpoly(ζ_X_,(0,-12,6,0,1)) = (ζ_X_^4+6*ζ_X_^2-12*ζ_X_)
    k_1 = evalpoly(ζ_X_,(0,-12,6,0,1))/(2*ζX3)
    k_2 = -3*ζ_X_^2/(8*(1-ζ_X_)^2)
    #(-ζ_X_^4+3*ζ_X_^2+3*ζ_X_) = evalpoly(ζ_X_,(0,3,3,0,-1))
    k_3 = evalpoly(ζ_X_,(0,3,3,0,-1))/(6*ζX3)
    return exp(k_0+x_0ij*k_1+x_0ij^2*k_2+x_0ij^3*k_3)
end

function ζeff_fdf(model::SAFTVRMieModel, V, T, z, λ,ζ_X_,ρ_S_)
    A = SAFTγMieconsts.A
    λ⁻¹ = one(λ)/λ
    Aλ⁻¹ = A * SA[one(λ); λ⁻¹; λ⁻¹*λ⁻¹; λ⁻¹*λ⁻¹*λ⁻¹]
    _f =  Aλ⁻¹ ⋅ SA[ζ_X_; ζ_X_^2; ζ_X_^3; ζ_X_^4]
    _df = Aλ⁻¹ ⋅  SA[1; 2ζ_X_; 3ζ_X_^2; 4ζ_X_^3] * ζ_X_/ρ_S_
    return _f,_df
end

function aS_1_fdf(model::SAFTVRMieModel, V, T, z, λ,ζ_X_= @f(ζ_X),ρ_S_ = @f(ρ_S))
    ζeff_,∂ζeff_ = @f(ζeff_fdf,λ,ζ_X_,ρ_S_)
    ζeff3 = (1-ζeff_)^3
    _f =  -1/(λ-3)*(1-ζeff_/2)/ζeff3
    _df = -1/(λ-3)*((1-ζeff_/2)/ζeff3
    + ρ_S_*((3*(1-ζeff_/2)*(1-ζeff_)^2
    - 0.5*ζeff3)/ζeff3^2)*∂ζeff_)
    return _f,_df
end

function B_fdf(model::SAFTVRMieModel, V, T, z, λ, x_0,ζ_X_= @f(ζ_X),ρ_S_ = @f(ρ_S))
    I = (1-x_0^(3-λ))/(λ-3)
    J = (1-(λ-3)*x_0^(4-λ)+(λ-4)*x_0^(3-λ))/((λ-3)*(λ-4))
    ζX2 = (1-ζ_X_)^2
    ζX3 = (1-ζ_X_)^3
    ζX6 = ζX3*ζX3
 
    _f = I*(1-ζ_X_/2)/ζX3-9*J*ζ_X_*(ζ_X_+1)/(2*ζX3)
    _df = (((1-ζ_X_/2)*I/ζX3-9*ζ_X_*(1+ζ_X_)*J/(2*ζX3))
        + ζ_X_*( (3*(1-ζ_X_/2)*ζX2
        - 0.5*ζX3)*I/ζX6
        - 9*J*((1+2*ζ_X_)*ζX3
        + ζ_X_*(1+ζ_X_)*3*ζX2)/(2*ζX6)))
    return _f,_df
end

function KHS_fdf(model::SAFTVRMieModel, V, T, z,ζ_X_,ρ_S_ = @f(ρ_S))
    ζX4 = (1-ζ_X_)^4
    denom1 = evalpoly(ζ_X_,(1,4,4,-4,1))
    ∂denom1 = evalpoly(ζ_X_,(4,8,-12,4))
    _f = ζX4/denom1
    _df = -(ζ_X_/ρ_S_)*((4*(1-ζ_X_)^3*denom1 + ζX4*∂denom1)/denom1^2)
    #@show _f,_df
    return _f,_df
end

function ∂a_2╱∂ρ_S(model::SAFTVRMieModel,V, T, z, i)
    λr = diagvalues(model.params.lambda_r)
    λa = diagvalues(model.params.lambda_a)
    x_0ij = @f(x_0,i,i)
    ζ_X_ = @f(ζ_X)
    ρ_S_ = @f(ρ_S)
    ∂KHS╱∂ρ_S = -ζ_X_/ρ_S_ *
    ( (4*(1-ζ_X_)^3*(1+4*ζ_X_+4*ζ_X_^2-4*ζ_X_^3+ζ_X_^4)
        + (1-ζ_X_)^4*(4+8*ζ_X_-12*ζ_X_^2+4*ζ_X_^3))/(1+4*ζ_X_+4*ζ_X_^2-4*ζ_X_^3+ζ_X_^4)^2 )
    return 0.5*@f(C,i,i)^2 *
    (@f(ρ_S)*∂KHS╱∂ρ_S*(x_0ij^(2*λa[i])*(@f(aS_1,2*λa[i])+@f(B,2*λa[i],x_0ij))
                         - 2*x_0ij^(λa[i]+λr[i])*(@f(aS_1,λa[i]+λr[i])+@f(B,λa[i]+λr[i],x_0ij))
                         + x_0ij^(2*λr[i])*(@f(aS_1,2*λr[i])+@f(B,2*λr[i],x_0ij)))
        + @f(KHS)*(x_0ij^(2*λa[i])*(@f(∂aS_1╱∂ρ_S,2*λa[i])+@f(∂B╱∂ρ_S,2*λa[i],x_0ij))
              - 2*x_0ij^(λa[i]+λr[i])*(@f(∂aS_1╱∂ρ_S,λa[i]+λr[i])+@f(∂B╱∂ρ_S,λa[i]+λr[i],x_0ij))
              + x_0ij^(2*λr[i])*(@f(∂aS_1╱∂ρ_S,2*λr[i])+@f(∂B╱∂ρ_S,2*λr[i],x_0ij))))
end

function I(model::SAFTVRMieModel, V, T, z,TR,_data = @f(data))
    _d,ρS,ζi,_ζ_X,_ζst,σ3_x = _data
    c  = SAFTVRMieconsts.c
    res = zero(_ζst)
    ρR = ρS*σ3_x
    @inbounds for n ∈ 0:10
        for m ∈ 0:(10-n)
            res += c[n+1,m+1]*ρR^n*TR^m
        end
    end
    return res
end

function Δ(model::SAFTVRMieModel, V, T, z, i, j, a, b,_data = @f(data))
    ϵ = model.params.epsilon
    TR = T/ϵ[i,j]
    _I = @f(I,TR,_data)
    ϵ_assoc = model.params.epsilon_assoc.values
    K = model.params.bondvol.values
    F = (exp(ϵ_assoc[i,j][a,b]/T)-1)
    return F*K[i,j][a,b]*_I
end

#optimized functions for maximum speed on default SAFTVRMie
function a_dispchain(model::SAFTVRMie, V, T, z,_data = @f(data))
    _d,ρS,ζi,_ζ_X,_ζst,_,m̄ = _data
    comps = @comps
    ∑z = ∑(z)
    m = model.params.segment
    _ϵ = model.params.epsilon
    _λr = model.params.lambda_r
    _λa = model.params.lambda_a
    _σ = model.params.sigma
    m̄inv = 1/m̄
    a₁ = zero(V+T+first(z))
    a₂ = a₁
    a₃ = a₁
    achain = a₁
    _ζst5 = _ζst^5
    _ζst8 = _ζst^8
    _KHS,_∂KHS = @f(KHS_fdf,_ζ_X,ρS)
    for i ∈ comps
        j = i
        x_Si = z[i]*m[i]*m̄inv
        x_Sj = x_Si
        ϵ = _ϵ[i,j]
        λa = _λa[i,j]
        λr = _λr[i,j] 
        σ = _σ[i,j]
        _C = @f(Cλ,λa,λr)
        dij = _d[i]
        x_0ij = σ/dij
        dij3 = dij^3
        #calculations for a1 - diagonal
        aS_1_a,∂aS_1∂ρS_a = @f(aS_1_fdf,λa,_ζ_X,ρS)
        aS_1_r,∂aS_1∂ρS_r = @f(aS_1_fdf,λr,_ζ_X,ρS)
        B_a,∂B∂ρS_a = @f(B_fdf,λa,x_0ij,_ζ_X,ρS)
        B_r,∂B∂ρS_r = @f(B_fdf,λr,x_0ij,_ζ_X,ρS)
        a1_ij = (2*π*ϵ*dij3)*_C*ρS*
        (x_0ij^λa*(aS_1_a+B_a) - x_0ij^λr*(aS_1_r+B_r))

        #calculations for a2 - diagonal
        aS_1_2a,∂aS_1∂ρS_2a = @f(aS_1_fdf,2*λa,_ζ_X,ρS)
        aS_1_2r,∂aS_1∂ρS_2r = @f(aS_1_fdf,2*λr,_ζ_X,ρS)
        aS_1_ar,∂aS_1∂ρS_ar = @f(aS_1_fdf,λa+λr,_ζ_X,ρS)
        B_2a,∂B∂ρS_2a = @f(B_fdf,2*λa,x_0ij,_ζ_X,ρS)
        B_2r,∂B∂ρS_2r = @f(B_fdf,2*λr,x_0ij,_ζ_X,ρS)
        B_ar,∂B∂ρS_ar = @f(B_fdf,λr+λa,x_0ij,_ζ_X,ρS)
        α = _C*(1/(λa-3)-1/(λr-3))
        f1,f2,f3,f4,f5,f6 = @f(f123456,α)
        _χ = f1*_ζst+f2*_ζst5+f3*_ζst8
        a2_ij = π*_KHS*(1+_χ)*ρS*ϵ^2*dij3*_C^2 *
        (x_0ij^(2*λa)*(aS_1_2a+B_2a)
        - 2*x_0ij^(λa+λr)*(aS_1_ar+B_ar)
        + x_0ij^(2*λr)*(aS_1_2r+B_2r))
        
        #calculations for a3 - diagonal
        a3_ij = -ϵ^3*f4*_ζst * exp(f5*_ζst+f6*_ζst^2)
        #adding - diagonal
        a₁ += a1_ij*x_Si*x_Sj
        a₂ += a2_ij*x_Si*x_Sj
        a₃ += a3_ij*x_Si*x_Sj

        g_HSi = @f(g_HS,x_0ij,_ζ_X)

        ∂a_1∂ρ_S = _C*(x_0ij^λa*(∂aS_1∂ρS_a+∂B∂ρS_a)
                      - x_0ij^λr*(∂aS_1∂ρS_r+∂B∂ρS_r))
        #calculus for g1
        g_1_ = 3*∂a_1∂ρ_S-_C*(λa*x_0ij^λa*(aS_1_a+B_a)-λr*x_0ij^λr*(aS_1_r+B_r))
        θ = exp(ϵ/T)-1
        γc =  10 * (-tanh(10*(0.57-α))+1) * _ζst*θ*exp(-6.7*_ζst-8*_ζst^2)
        ∂a_2∂ρ_S = 0.5*_C^2 *
            (ρS*_∂KHS*(x_0ij^(2*λa)*(aS_1_2a+B_2a)
            - 2*x_0ij^(λa+λr)*(aS_1_ar+B_ar)
            + x_0ij^(2*λr)*(aS_1_2r+B_2r))
            + _KHS*(x_0ij^(2*λa)*(∂aS_1∂ρS_2a+∂B∂ρS_2a)
            - 2*x_0ij^(λa+λr)*(∂aS_1∂ρS_ar+∂B∂ρS_ar)
            + x_0ij^(2*λr)*(∂aS_1∂ρS_2r+∂B∂ρS_2r)))
    
        gMCA2 = 3*∂a_2∂ρ_S-_KHS*_C^2 *
        (λr*x_0ij^(2*λr)*(aS_1_2r+B_2r)-
            (λa+λr)*x_0ij^(λa+λr)*(aS_1_ar+B_ar)+
            λa*x_0ij^(2*λa)*(aS_1_2a+B_2a))
        g_2_ = (1+γc)*gMCA2
        g_Mie_ = g_HSi*exp(ϵ/T*g_1_/g_HSi+(ϵ/T)^2*g_2_/g_HSi)
        achain -=  z[i]*(log(g_Mie_)*(m[i]-1))
        for j ∈ 1:i-1
            x_Sj = z[j]*m[j]*m̄inv   
            ϵ = _ϵ[i,j]
            λa = _λa[i,j]
            λr = _λr[i,j] 
            σ = _σ[i,j]
            _C = @f(Cλ,λa,λr)
            dij = 0.5*(_d[i]+_d[j])
            x_0ij = σ/dij
            dij3 = dij^3
            #calculations for a1
            a1_ij = (2*π*ϵ*dij3)*_C*ρS*
            (x_0ij^λa*(@f(aS_1,λa,_ζ_X)+@f(B,λa,x_0ij,_ζ_X)) - x_0ij^λr*(@f(aS_1,λr,_ζ_X)+@f(B,λr,x_0ij,_ζ_X)))

            #calculations for a2
            α = _C*(1/(λa-3)-1/(λr-3))
            f1,f2,f3,f4,f5,f6 = @f(f123456,α)
            _χ = f1*_ζst+f2*_ζst5+f3*_ζst8
            a2_ij = π*_KHS*(1+_χ)*ρS*ϵ^2*dij3*_C^2 *
            (x_0ij^(2*λa)*(@f(aS_1,2*λa,_ζ_X)+@f(B,2*λa,x_0ij,_ζ_X))
            - 2*x_0ij^(λa+λr)*(@f(aS_1,λa+λr,_ζ_X)+@f(B,λa+λr,x_0ij,_ζ_X))
            + x_0ij^(2*λr)*(@f(aS_1,2λr,_ζ_X)+@f(B,2*λr,x_0ij,_ζ_X)))
            
            #calculations for a3
            a3_ij = -ϵ^3*f4*_ζst * exp(f5*_ζst+f6*_ζst^2)
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
    return adisp + achain/∑z
end

function a_disp(model::SAFTVRMieModel, V, T, z,_data = @f(data))
    _d,ρS,ζi,_ζ_X,_ζst,_,m̄ = _data 
    comps = 1:length(z)
    #this is a magic trick. we normally (should) expect length(z) = length(model),
    #but on GC models, @comps != @groups
    #if we pass Xgc instead of z, the equation is exactly the same.
    #we need to add the divide the result by sum(z) later.
    m = model.params.segment
    _ϵ = model.params.epsilon
    _λr = model.params.lambda_r
    _λa = model.params.lambda_a
    _σ = model.params.sigma
    m̄inv = 1/m̄
    a₁ = zero(V+T+first(z))
    a₂ = a₁
    a₃ = a₁
    _ζst5 = _ζst^5
    _ζst8 = _ζst^8
    _KHS = @f(KHS,_ζ_X,ρS)
    for i ∈ comps
        j = i
        x_Si = z[i]*m[i]*m̄inv
        x_Sj = x_Si
        ϵ = _ϵ[i,j]
        λa = _λa[i,i]
        λr = _λr[i,i] 
        σ = _σ[i,i]
        _C = @f(Cλ,λa,λr)
        dij = _d[i]
        dij3 = dij^3
        x_0ij = σ/dij
        #calculations for a1 - diagonal
        aS_1_a = @f(aS_1,λa,_ζ_X)
        aS_1_r = @f(aS_1,λr,_ζ_X)
        B_a = @f(B,λa,x_0ij,_ζ_X)
        B_r = @f(B,λr,x_0ij,_ζ_X)
        a1_ij = (2*π*ϵ*dij3)*_C*ρS*
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
        a2_ij = π*_KHS*(1+_χ)*ρS*ϵ^2*dij3*_C^2 *
        (x_0ij^(2*λa)*(aS_1_2a+B_2a)
        - 2*x_0ij^(λa+λr)*(aS_1_ar+B_ar)
        + x_0ij^(2*λr)*(aS_1_2r+B_2r))
        
        #calculations for a3 - diagonal
        a3_ij = -ϵ^3*f4*_ζst * exp(f5*_ζst+f6*_ζst^2)
        #adding - diagonal
        a₁ += a1_ij*x_Si*x_Si
        a₂ += a2_ij*x_Si*x_Si
        a₃ += a3_ij*x_Si*x_Si
        for j ∈ 1:(i-1)
            x_Sj = z[j]*m[j]*m̄inv   
            ϵ = _ϵ[i,j]
            λa = _λa[i,j]
            λr = _λr[i,j] 
            σ = _σ[i,j]
            _C = @f(Cλ,λa,λr)
            dij = 0.5*(_d[i]+_d[j])
            x_0ij = σ/dij
            dij3 = dij^3
            x_0ij = σ/dij
            #calculations for a1
            a1_ij = (2*π*ϵ*dij3)*_C*ρS*
            (x_0ij^λa*(@f(aS_1,λa,_ζ_X)+@f(B,λa,x_0ij,_ζ_X)) - x_0ij^λr*(@f(aS_1,λr,_ζ_X)+@f(B,λr,x_0ij,_ζ_X)))

            #calculations for a2
            α = _C*(1/(λa-3)-1/(λr-3))
            f1,f2,f3,f4,f5,f6 = @f(f123456,α)
            _χ = f1*_ζst+f2*_ζst5+f3*_ζst8
            a2_ij = π*_KHS*(1+_χ)*ρS*ϵ^2*dij3*_C^2 *
            (x_0ij^(2*λa)*(@f(aS_1,2*λa,_ζ_X)+@f(B,2*λa,x_0ij,_ζ_X))
            - 2*x_0ij^(λa+λr)*(@f(aS_1,λa+λr,_ζ_X)+@f(B,λa+λr,x_0ij,_ζ_X))
            + x_0ij^(2*λr)*(@f(aS_1,2λr,_ζ_X)+@f(B,2*λr,x_0ij,_ζ_X)))
            
            #calculations for a3
            a3_ij = -ϵ^3*f4*_ζst * exp(f5*_ζst+f6*_ζst^2)
            #adding
            a₁ += 2*a1_ij*x_Si*x_Sj
            a₂ += 2*a2_ij*x_Si*x_Sj
            a₃ += 2*a3_ij*x_Si*x_Sj            
        end
    end
    a₁ = a₁*m̄/T #/sum(z)
    a₂ = a₂*m̄/(T*T)  #/sum(z)
    a₃ = a₃*m̄/(T*T*T)  #/sum(z)
    #@show (a₁,a₂,a₃)
    adisp =  a₁ + a₂ + a₃ 
    return adisp
end

function a_chain(model::SAFTVRMieModel, V, T, z,_data = @f(data))
    _d,ρS,ζi,_ζ_X,_ζst,_,m̄ = _data
    l = length(z)
    comps = 1:l
    ∑z = ∑(z)
    m = model.params.segment
    _ϵ = model.params.epsilon
    _λr = model.params.lambda_r
    _λa = model.params.lambda_a
    _σ = model.params.sigma
    m̄inv = 1/m̄
    a₁ = zero(V+T+first(z))
    a₂ = a₁
    a₃ = a₁
    achain = a₁
    _ζst5 = _ζst^5
    _ζst8 = _ζst^8
    _KHS,_∂KHS = @f(KHS_fdf,_ζ_X,ρS)
    for i ∈ comps
        x_Si = z[i]*m[i]*m̄inv
        x_Sj = x_Si
        ϵ = _ϵ[i,i]
        λa = _λa[i,i]
        λr = _λr[i,i] 
        σ = _σ[i,i]
        _C = @f(Cλ,λa,λr)
        dij = _d[i]
        x_0ij = σ/dij
        dij3 = dij^3
        x_0ij = σ/dij
        #calculations for a1 - diagonal
        aS_1_a,∂aS_1∂ρS_a = @f(aS_1_fdf,λa,_ζ_X,ρS)
        aS_1_r,∂aS_1∂ρS_r = @f(aS_1_fdf,λr,_ζ_X,ρS)
        B_a,∂B∂ρS_a = @f(B_fdf,λa,x_0ij,_ζ_X,ρS)
        B_r,∂B∂ρS_r = @f(B_fdf,λr,x_0ij,_ζ_X,ρS)
        a1_ij = (2*π*ϵ*dij3)*_C*ρS*
        (x_0ij^λa*(aS_1_a+B_a) - x_0ij^λr*(aS_1_r+B_r))

        #calculations for a2 - diagonal
        aS_1_2a,∂aS_1∂ρS_2a = @f(aS_1_fdf,2*λa,_ζ_X,ρS)
        aS_1_2r,∂aS_1∂ρS_2r = @f(aS_1_fdf,2*λr,_ζ_X,ρS)
        aS_1_ar,∂aS_1∂ρS_ar = @f(aS_1_fdf,λa+λr,_ζ_X,ρS)
        B_2a,∂B∂ρS_2a = @f(B_fdf,2*λa,x_0ij,_ζ_X,ρS)
        B_2r,∂B∂ρS_2r = @f(B_fdf,2*λr,x_0ij,_ζ_X,ρS)
        B_ar,∂B∂ρS_ar = @f(B_fdf,λr+λa,x_0ij,_ζ_X,ρS)
        α = _C*(1/(λa-3)-1/(λr-3))
        f1,f2,f3,f4,f5,f6 = @f(f123456,α)
        _χ = f1*_ζst+f2*_ζst5+f3*_ζst8
        a2_ij = π*_KHS*(1+_χ)*ρS*ϵ^2*dij3*_C^2 *
        (x_0ij^(2*λa)*(aS_1_2a+B_2a)
        - 2*x_0ij^(λa+λr)*(aS_1_ar+B_ar)
        + x_0ij^(2*λr)*(aS_1_2r+B_2r))
        
        #calculations for a3 - diagonal
        a3_ij = -ϵ^3*f4*_ζst * exp(f5*_ζst+f6*_ζst^2)
        #adding - diagonal
        a₁ += a1_ij*x_Si*x_Sj
        a₂ += a2_ij*x_Si*x_Sj
        a₃ += a3_ij*x_Si*x_Sj

        g_HSi = @f(g_HS,x_0ij,_ζ_X)
        #@show (g_HSi,i)
        ∂a_1∂ρ_S = _C*(x_0ij^λa*(∂aS_1∂ρS_a+∂B∂ρS_a)
                      - x_0ij^λr*(∂aS_1∂ρS_r+∂B∂ρS_r))
        #@show (∂a_1∂ρ_S,1)

        g_1_ = 3*∂a_1∂ρ_S-_C*(λa*x_0ij^λa*(aS_1_a+B_a)-λr*x_0ij^λr*(aS_1_r+B_r))
        #@show (g_1_,i)
        θ = exp(ϵ/T)-1
        γc =  10 * (-tanh(10*(0.57-α))+1) * _ζst*θ*exp(-6.7*_ζst-8*_ζst^2)
        ∂a_2∂ρ_S = 0.5*_C^2 *
            (ρS*_∂KHS*(x_0ij^(2*λa)*(aS_1_2a+B_2a)
            - 2*x_0ij^(λa+λr)*(aS_1_ar+B_ar)
            + x_0ij^(2*λr)*(aS_1_2r+B_2r))
            + _KHS*(x_0ij^(2*λa)*(∂aS_1∂ρS_2a+∂B∂ρS_2a)
            - 2*x_0ij^(λa+λr)*(∂aS_1∂ρS_ar+∂B∂ρS_ar)
            + x_0ij^(2*λr)*(∂aS_1∂ρS_2r+∂B∂ρS_2r)))
    
        gMCA2 = 3*∂a_2∂ρ_S-_KHS*_C^2 *
        (λr*x_0ij^(2*λr)*(aS_1_2r+B_2r)-
            (λa+λr)*x_0ij^(λa+λr)*(aS_1_ar+B_ar)+
            λa*x_0ij^(2*λa)*(aS_1_2a+B_2a))
        g_2_ = (1+γc)*gMCA2
        #@show (g_2_,i)
        g_Mie_ = g_HSi*exp(ϵ/T*g_1_/g_HSi+(ϵ/T)^2*g_2_/g_HSi)
        #@show (g_Mie_,i)
        achain +=  z[i]*(log(g_Mie_)*(m[i]-1))
    end
    return -achain/∑z
end
const SAFTVRMieconsts = (
    u = (0.26356031971814109102031,1.41340305910651679221800,3.59642577104072208122300,7.08581000585883755692200,12.6408008442757826594300),
    w = (0.5217556105828086524759,0.3986668110831759274500,7.5942449681707595390e-2,3.6117586799220484545e-3,2.3369972385776227891e-5),

    A = SA[0.81096   1.7888  -37.578   92.284;
    1.02050  -19.341   151.26  -463.50;
    -1.90570   22.845  -228.14   973.92;
    1.08850  -6.1962   106.98  -677.64],

    ϕ = [(7.5365557, -359.440,  1550.9, -1.199320, -1911.2800,  9236.9),
        (-37.604630,  1825.60, -5070.1,  9.063632,  21390.175, -129430.0),
        (71.745953, -3168.00,  6534.6, -17.94820, -51320.700,  357230.0),
        (-46.835520,  1884.20, -3288.7,  11.34027,  37064.540, -315530.0),
        (-2.4679820,- 0.82376, -2.7171,  20.52142,  1103.7420,  1390.2),
        (-0.5027200, -3.19350,  2.0883, -56.63770, -3264.6100, -4518.2),
        (8.0956883,  3.70900,  0.0000,  40.53683,  2556.1810,  4241.6)],

    c  = [0.0756425183020431	-0.128667137050961	 0.128350632316055	-0.0725321780970292	   0.0257782547511452  -0.00601170055221687	  0.000933363147191978  -9.55607377143667e-05  6.19576039900837e-06 -2.30466608213628e-07 3.74605718435540e-09
          0.134228218276565	    -0.182682168504886 	 0.0771662412959262	-0.000717458641164565 -0.00872427344283170	0.00297971836051287	 -0.000484863997651451	 4.35262491516424e-05 -2.07789181640066e-06	4.13749349344802e-08 0
         -0.565116428942893	     1.00930692226792   -0.660166945915607	 0.214492212294301	  -0.0388462990166792	0.00406016982985030	 -0.000239515566373142	 7.25488368831468e-06 -8.58904640281928e-08	0	                 0
         -0.387336382687019	    -0.211614570109503	 0.450442894490509	-0.176931752538907	   0.0317171522104923  -0.00291368915845693	  0.000130193710011706  -2.14505500786531e-06  0	                0	                 0
          2.13713180911797	    -2.02798460133021 	 0.336709255682693	 0.00118106507393722  -0.00600058423301506	0.000626343952584415 -2.03636395699819e-05	 0	                   0	                0	                 0
         -0.300527494795524	     2.89920714512243   -0.567134839686498	 0.0518085125423494	  -0.00239326776760414	4.15107362643844e-05  0	                     0	                   0	                0                    0
         -6.21028065719194	    -1.92883360342573	 0.284109761066570	-0.0157606767372364	   0.000368599073256615	0 	                  0	                     0	                   0	                0	                 0
          11.6083532818029	     0.742215544511197  -0.0823976531246117	 0.00186167650098254   0	                0	                  0	                     0	                   0	                0	                 0
         -10.2632535542427	    -0.125035689035085	 0.0114299144831867	 0	                   0	                0	                  0	                     0	                   0	                0	                 0
          4.65297446837297	    -0.00192518067137033 0	                 0	                   0	                0	                  0	                     0	                   0	                0	                 0
         -0.867296219639940	     0	                 0	                 0	                   0	                0	                  0	                     0	                   0	                0	                 0],
)

########
#=
Optimizations for single component SAFTVRMie
=#

#######

function d(model::SAFTVRMie, V, T, z::SingleComp)
    _ϵ = model.params.epsilon
    _σ = model.params.sigma
    _λa = model.params.lambda_a
    _λr = model.params.lambda_r
    u = SAFTVRMieconsts.u
    w = SAFTVRMieconsts.w
    ϵ = _ϵ[1]
    σ = _σ[1]
    λa = _λa[1]
    λr = _λr[1]
    θ = Cλ(model,V,T,z,λa,λr)*ϵ/T
    λrinv = 1/λr
    λaλr = λa/λr
    di = zero(T*1.0)
    for j ∈ 1:5
        di += w[j]*(θ/(θ+u[j]))^λrinv*(exp(θ*(1/(θ/(θ+u[j]))^λaλr-1))/(u[j]+θ)/λr)
    end
    return SA[σ*(1-di)]
end
