struct SAFTVRMieParam{T} <: ParametricEoSParam{T}
    Mw::SingleParam{T}
    segment::SingleParam{T}
    sigma::PairParam{T}
    lambda_a::PairParam{T}
    lambda_r::PairParam{T}
    epsilon::PairParam{T}
    epsilon_assoc::AssocParam{T}
    bondvol::AssocParam{T}
end

function SAFTVRMieParam(Mw,segment,sigma,lambda_a,lambda_r,epsilon,epsilon_assoc,bondvol)
    return build_parametric_param(SAFTVRMieParam,Mw,segment,sigma,lambda_a,lambda_r,epsilon,epsilon_assoc,bondvol) 
end

abstract type SAFTVRMieModel <: SAFTModel end
@newmodel SAFTVRMie SAFTVRMieModel SAFTVRMieParam{T}
default_references(::Type{SAFTVRMie}) = ["10.1063/1.4819786", "10.1080/00268976.2015.1029027"]
default_locations(::Type{SAFTVRMie}) = ["SAFT/SAFTVRMie", "properties/molarmass.csv"]
function transform_params(::Type{SAFTVRMie},params)
    sigma = params["sigma"]
    sigma.values .*= 1E-10
    sigma = sigma_LorentzBerthelot(sigma)
    epsilon = epsilon_HudsenMcCoubreysqrt(params["epsilon"], sigma)
    lambda_a = lambda_LorentzBerthelot(params["lambda_a"])
    lambda_r = lambda_LorentzBerthelot(params["lambda_r"])
    params["sigma"] = sigma
    params["epsilon"] = epsilon
    params["lambda_a"] = lambda_a
    params["lambda_r"] = lambda_r
    return params
end
"""
    SAFTVRMieModel <: SAFTModel

    SAFTVRMie(components;
    idealmodel = BasicIdeal,
    userlocations = String[],
    ideal_userlocations = String[],
    reference_state = nothing,
    verbose = false,
    assoc_options = AssocOptions())

## Input parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Single Parameter (`Float64`) - Segment Diameter `[Å]`
- `epsilon`: Single Parameter (`Float64`) - Reduced dispersion energy `[K]`
- `lambda_a`: Pair Parameter (`Float64`) - Atractive range parameter (no units)
- `lambda_r`: Pair Parameter (`Float64`) - Repulsive range parameter (no units)
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Parameter (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m³]`

## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Pair Parameter (`Float64`) - Mixed segment Diameter `[m]`
- `lambda_a`: Pair Parameter (`Float64`) - Atractive range parameter (no units)
- `lambda_r`: Pair Parameter (`Float64`) - Repulsive range parameter (no units)
- `epsilon`: Pair Parameter (`Float64`) - Mixed reduced dispersion energy`[K]`
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m³]`

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

include("equations.jl")

function recombine_impl!(model::SAFTVRMieModel)
    sigma = model.params.sigma
    epsilon = model.params.epsilon
    lambda_a = model.params.lambda_a
    lambda_r = model.params.lambda_r
    sigma = sigma_LorentzBerthelot!(sigma)
    epsilon = epsilon_HudsenMcCoubreysqrt!(epsilon,sigma)
    lambda_a = lambda_LorentzBerthelot!(lambda_a)
    lambda_r = lambda_LorentzBerthelot!(lambda_r)
    recombine_assoc!(model)
    return model
end

function x0_volume_liquid(model::SAFTVRMieModel,T,z)
    v_lb = lb_volume(model,T,z)
    return v_lb*1.5
end

function data(model::SAFTVRMieModel, V, T, z)
    m̄ = dot(z,model.params.segment.values)
    _d = @f(d)
    ζi = @f(ζ0123,_d)
    _ζ_X,σ3x = @f(ζ_X_σ3,_d,m̄)
    _ρ_S = @f(ρ_S,m̄)
    _ζst = σ3x*_ρ_S*π/6
    return (_d,_ρ_S,ζi,_ζ_X,_ζst,σ3x,m̄)
end

function packing_fraction(model::SAFTVRMieModel,_data)
    _,_,ζi,_,_,_,m̄ = _data
    _,_,_,η = ζi
    return η
end

# function a_res(model::SAFTVRMieModel, V, T, z, _data = @f(data))
#     return @f(a_hs,_data)+@f(a_disp,_data) + @f(a_chain,_data) + @f(a_assoc,_data)
# end

#fused chain and disp calculation
function a_res(model::SAFTVRMieModel, V, T, z, _data = @f(data))
    return @f(a_hs,_data) + @f(a_dispchain,_data) + @f(a_assoc,_data)
end

function a_mono(model::SAFTVRMieModel, V, T, z,_data = @f(data))
    return @f(a_hs,_data)+@f(a_disp,_data)
end

function a_hs(model::SAFTVRMieModel, V, T, z,_data = @f(data))
    _d,_,ζi,_,_,_,m̄ = _data 
    ζ0,ζ1,ζ2,ζ3 = ζi
    if !iszero(ζ3)
        _a_hs = bmcs_hs(ζ0,ζ1,ζ2,ζ3)
    else
        _a_hs = @f(bmcs_hs_zero_v,_d)
    end

    return m̄*_a_hs/sum(z)
end

function ρ_S(model::SAFTVRMieModel, V, T, z, m̄ = dot(z,model.params.segment.values))
    return N_A/V*m̄
end

#=
SAFT-VR-Mie diameter:
Defined as:
```
C  = (λr/(λr-λa))*(λr/λa)^(λa/(λr-λa))
u(r) = C*ϵ*(x^-λr - x^-λa)
f(r) = exp(-u(r)/T)
d = σ*(1-integral(f(r),0,1))
```

we use a mixed approach, depending on T⋆ = T/ϵ:

if T⋆ < 1:
    5-point gauss-laguerre. we do the change of variables `y = r^-λr`
else:
    10-point modified gauss-legendre with cut.
=#
function d_vrmie(T,λa,λr,σ,ϵ)
    C = Cλ_mie(λa, λr)
    θ = ϵ*C/T
    ∑fi = vr_mie_d_integral(θ,λa,λr)
    return σ*(1 - ∑fi)
end

#this function is a fundamental one.
function vr_mie_d_integral(θ,λa,λr)
    λrinv = 1/λr
    λaλr = λa/λr
    if θ > 1
        function f_laguerre(x)
            lnx = log(x)
            return exp(-λrinv*lnx)*exp(θ*exp(lnx*λaλr))*λrinv/x
        end
        return Solvers.laguerre10(f_laguerre,θ,one(θ))
    else
        j = d_vrmie_cut(θ,λa,λr)
        function f_legendre(x) 
            lnx = log(x)
            return exp(-θ*(exp(-λr*lnx)-exp(-λa*lnx)))
        end
        return Solvers.integral10(f_legendre,j,one(j))
    end
end

#implements the method of aasen for VRQ Mie. (https://github.com/usnistgov/teqp/issues/39)
function d_vrmie_cut(θ,λa,λr)
    #initial point
    EPS = eps(typeof(θ))
    K = log(-log(EPS)/θ)
    j0 = exp(-K/λr)
    # exp(-u(r)/T), d[exp(-u(r))/T)]/dr, d2[exp(-u(r))/T)]/dr2
    function fdfd2f(r)
        r⁻¹ = 1/r
        lnr = log(r)
        rλr = exp(-lnr*λr)
        rλa = exp(-lnr*λa)
        #rλr = r^-λr
        #rλa = r^-λa
        u_r = rλr - rλa #u/C*ϵ
        du_ra = rλa*r⁻¹*(-λa)
        du_rr = rλr*r⁻¹*(-λr)
        du_r = (du_rr - du_ra)
        d2u_rr = du_rr*r⁻¹*(-λr - 1)
        d2u_ra = du_ra*r⁻¹*(-λa - 1)
        d2u_r = d2u_rr - d2u_ra
        f = exp(-u_r*θ)
        df = -θ*f*du_r
        d2f = df*df - θ*d2u_r*f
        return f, f/df, df/d2f
    end
    j = j0
    for i in 1:5
        fi,f1,f2 = fdfd2f(j)
        dd = (1 - 0.5*f1/f2)
        dj = f1/(1 - 0.5*f1/f2)
        j = j - dj
        fi < eps(eltype(fi)) && break 
    end
    return j
end

function d(model::SAFTVRMieModel, V, T, z)
    ϵ = diagvalues(model.params.epsilon.values)
    σ = diagvalues(model.params.sigma.values)
    λa = diagvalues(model.params.lambda_a.values)
    λr = diagvalues(model.params.lambda_r.values)
    n = length(z)
    _d = fill(zero(V+T+first(z)+one(eltype(model))),n)
    for k ∈ 1:n
        _d[k] = d_vrmie(T,λa[k],λr[k],σ[k],ϵ[k])
    end
    return _d
end


function d(model::SAFTVRMieModel, V, T, z, λa,λr,ϵ,σ)
    d_vrmie(T,λa,λr,σ,ϵ)
end

function Cλ(model::SAFTVRMieModel, V, T, z, λa, λr)
    return Cλ_mie(λa, λr)
end

Cλ_mie(λa, λr) = (λr/(λr-λa))*(λr/λa)^(λa/(λr-λa))

function ζ_X(model::SAFTVRMieModel, V, T, z,_d = @f(d))
    _ζ_X,σ3x = @f(ζ_X_σ3,_d)
    return _ζ_X
end

function ζ_X_σ3(model::SAFTVRMieModel, V, T, z,_d = @f(d),m̄ = dot(z,model.params.segment.values))
    m = model.params.segment.values
    m̄ = dot(z, m)
    m̄inv = 1/m̄
    σ = model.params.sigma.values
    ρS = N_A/V*m̄
    comps = 1:length(z)
    _ζ_X = zero(V+T+first(z)+one(eltype(model)))
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

aS_1(model::SAFTVRMieModel, V, T, z, λ,ζ_X_= @f(ζ_X)) = SAFTVRMieUtils.aS1(λ,ζ_X_)

ζeff(model::SAFTVRMieModel, V, T, z, λ,ζ_X_= @f(ζ_X)) = SAFTVRMieUtils.ζeff(λ,ζ_X_)

B(model::SAFTVRMieModel, V, T, z, λ, x_0,ζ_X_ = @f(ζ_X)) = SAFTVRMieUtils.B(λ, x_0, ζ_X_)

KHS(model::SAFTVRMieModel, V, T, z,ζ_X_ = @f(ζ_X),ρS = nothing) = SAFTVRMieUtils.KHS(ζ_X_)

f123456(model::SAFTVRMieModel, V, T, z, α) = SAFTVRMieUtils.f123456(α)

function ζst(model::SAFTVRMieModel, V, T, z, _σ = model.params.sigma.values)
    m = model.params.segment.values
    m̄ = dot(z, m)
    m̄inv = 1/m̄
    ρS = N_A/V*m̄
    comps = @comps
    _ζst = zero(V+T+first(z)+one(eltype(model)))
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

g_HS(model::SAFTVRMieModel, V, T, z, x_0ij,ζ_X_ = @f(ζ_X)) = SAFTVRMieUtils.g_HS(x_0ij,ζ_X_)

function ζeff_fdf(model::SAFTVRMieModel, V, T, z, λ,ζ_X_,ρ_S_)
    A = SAFTγMieconsts.A
    λ⁻¹ = one(λ)/λ
    Aλ⁻¹ = A * SA[one(λ); λ⁻¹; λ⁻¹*λ⁻¹; λ⁻¹*λ⁻¹*λ⁻¹]
    _f = dot(Aλ⁻¹,SA[ζ_X_; ζ_X_^2; ζ_X_^3; ζ_X_^4])
    _df = dot(Aλ⁻¹,SA[1; 2ζ_X_; 3ζ_X_^2; 4ζ_X_^3]) * ζ_X_/ρ_S_
    return _f,_df
end

function ζeff_f_ρdf(model::SAFTVRMieModel, V, T, z, λ,ζ_X_)
    A = SAFTγMieconsts.A
    λ⁻¹ = one(λ)/λ
    Aλ⁻¹ = A * SA[one(λ); λ⁻¹; λ⁻¹*λ⁻¹; λ⁻¹*λ⁻¹*λ⁻¹]
    _f = dot(Aλ⁻¹,SA[ζ_X_; ζ_X_^2; ζ_X_^3; ζ_X_^4])
    _ρdf = dot(Aλ⁻¹,SA[1; 2ζ_X_; 3ζ_X_^2; 4ζ_X_^3]) * ζ_X_
    return _f,_ρdf
end

function aS_1_fdf(model::SAFTVRMieModel, V, T, z, λ, ζ_X_= @f(ζ_X),ρ_S_ = 0.0)
    ζeff_,∂ζeff_ρ_S = @f(ζeff_f_ρdf,λ,ζ_X_)
    ζeff3 = (1-ζeff_)^3
    ζeffm1 = (1-ζeff_*0.5)
    ζf = ζeffm1/ζeff3
    λf = -1/(λ-3)
    _f = λf * ζf
    _df = λf * (ζf + ∂ζeff_ρ_S*((3*ζeffm1*(1-ζeff_)^2 - 0.5*ζeff3)/ζeff3^2))
    return _f,_df
end

function B_fdf(model::SAFTVRMieModel, V, T, z, λ, x_0,ζ_X_= @f(ζ_X),ρ_S_ = @f(ρ_S))
    x_0_λ = x_0^(3-λ)
    I = (1-x_0_λ)/(λ-3)
    J = (1-(λ-3)*x_0^(4-λ)+(λ-4)*x_0_λ)/((λ-3)*(λ-4))
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
    _f,_ρdf = KHS_f_ρdf(model,V,T,z,ζ_X_)
    _df = _ρdf/ρ_S_
    return _f,_ρdf/ρ_S_
end

function KHS_f_ρdf(model::SAFTVRMieModel, V, T, z,ζ_X_)
    ζX4 = (1-ζ_X_)^4
    denom1 = evalpoly(ζ_X_,(1,4,4,-4,1))
    ∂denom1 = evalpoly(ζ_X_,(4,8,-12,4))
    _f = ζX4/denom1
    _df = -ζ_X_*((4*(1-ζ_X_)^3*denom1 + ζX4*∂denom1)/denom1^2)

    return _f,_df
end

function ∂a_2╱∂ρ_S(model::SAFTVRMieModel,V, T, z, i)
    λr = diagvalues(model.params.lambda_r.values)
    λa = diagvalues(model.params.lambda_a.values)
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

function I(model::SAFTVRMieModel, V, T, z, i, j, _data = @f(data))
    ϵ = model.params.epsilon.values[i,j]
    Tr = T/ϵ
    _d,ρS,ζi,_ζ_X,_ζst,σ3_x = _data
    c  = SAFTVRMieconsts.c
    res = zero(_ζst)
    ρr = ρS*σ3_x
    ρrn = one(ρr)
    @inbounds for n ∈ 0:10
        res_m = zero(res)
        Trm = one(Tr)
        for m ∈ 0:(10-n)
            res_m += c[n+1,m+1]*Trm
            Trm = Trm*Tr
        end
        res += res_m*ρrn
        ρrn = ρrn*ρr
    end
    return res
end

function Δ(model::SAFTVRMieModel, V, T, z, i, j, a, b,_data = @f(data))
    K = model.params.bondvol.values
    Kijab = K[i,j][a,b]
    if iszero(Kijab)
        return zero(@f(Base.promote_eltype))
    end
    _I = @f(I,i,j,_data)
    ϵ_assoc = model.params.epsilon_assoc.values
    F = expm1(ϵ_assoc[i,j][a,b]/T)
    return F*Kijab*_I
end

#optimized functions for maximum speed on default SAFTVRMie
function a_dispchain(model::SAFTVRMieModel, V, T, z,_data = @f(data))
    _d,ρS,ζi,ζₓ,_ζst,_,m̄ = _data
    comps = @comps
    ∑z = ∑(z)
    m = model.params.segment.values
    _ϵ = model.params.epsilon.values
    _λr = model.params.lambda_r.values
    _λa = model.params.lambda_a.values
    _σ = model.params.sigma.values
    m̄inv = 1/m̄
    a₁ = zero(V+T+first(z)+one(eltype(model)))
    a₂ = a₁
    a₃ = a₁
    achain = a₁
    _ζst5 = _ζst^5
    _ζst8 = _ζst^8
    _KHS,ρS_∂KHS = @f(KHS_f_ρdf,ζₓ)
    for i ∈ comps
        j = i
        mi = m[i]
        x_Si = z[i]*mi*m̄inv
        x_Sj = x_Si
        ϵ = _ϵ[i,j]
        λa = _λa[i,j]
        λr = _λr[i,j]
        σ = _σ[i,j]
        _C = @f(Cλ,λa,λr)
        dij = _d[i]
        x_0ij = σ/dij
        dij3 = dij^3
        τ = ϵ/T
        #precalculate exponentials of x_0ij
        x_0ij_λa = x_0ij^λa
        x_0ij_λr = x_0ij^λr
        x_0ij_2λa = x_0ij^(2*λa)
        x_0ij_2λr = x_0ij^(2*λr)
        x_0ij_λaλr = x_0ij^(λa + λr)

        #calculations for a1 - diagonal
        aS₁_a,∂aS₁∂ρS_a = @f(aS_1_fdf,λa,ζₓ,ρS)
        aS₁_r,∂aS₁∂ρS_r = @f(aS_1_fdf,λr,ζₓ,ρS)
        B_a,∂B∂ρS_a = @f(B_fdf,λa,x_0ij,ζₓ,ρS)
        B_r,∂B∂ρS_r = @f(B_fdf,λr,x_0ij,ζₓ,ρS)
        a1_ij = (2*π*ϵ*dij3)*_C*ρS*
        (x_0ij_λa*(aS₁_a+B_a) - x_0ij_λr*(aS₁_r+B_r))

        #calculations for a2 - diagonal
        aS₁_2a,∂aS₁∂ρS_2a = @f(aS_1_fdf,2*λa,ζₓ,ρS)
        aS₁_2r,∂aS₁∂ρS_2r = @f(aS_1_fdf,2*λr,ζₓ,ρS)
        aS₁_ar,∂aS₁∂ρS_ar = @f(aS_1_fdf,λa+λr,ζₓ,ρS)
        B_2a,∂B∂ρS_2a = @f(B_fdf,2*λa,x_0ij,ζₓ,ρS)
        B_2r,∂B∂ρS_2r = @f(B_fdf,2*λr,x_0ij,ζₓ,ρS)
        B_ar,∂B∂ρS_ar = @f(B_fdf,λr+λa,x_0ij,ζₓ,ρS)
        α = _C*(1/(λa-3)-1/(λr-3))
        f1,f2,f3,f4,f5,f6 = @f(f123456,α)
        _χ = f1*_ζst + f2*_ζst5 + f3*_ζst8
        a2_ij = π*_KHS*(1+_χ)*ρS*ϵ^2*dij3*_C^2 *
                (x_0ij_2λa*(aS₁_2a+B_2a)
                    - 2*x_0ij_λaλr*(aS₁_ar+B_ar)
                    + x_0ij_2λr*(aS₁_2r+B_2r)
                )

        #calculations for a3 - diagonal
        a3_ij = -ϵ^3*f4*_ζst*exp(_ζst*(f5 + f6*_ζst))
        #adding - diagonal
        a₁ += a1_ij*x_Si*x_Sj
        a₂ += a2_ij*x_Si*x_Sj
        a₃ += a3_ij*x_Si*x_Sj

        g_HSi = @f(g_HS,x_0ij,ζₓ)

        ∂a_1∂ρ_S = _C*(x_0ij_λa*(∂aS₁∂ρS_a+∂B∂ρS_a)
                    - x_0ij_λr*(∂aS₁∂ρS_r+∂B∂ρS_r)
                    )
        #calculus for g1
        g_1_ = 3*∂a_1∂ρ_S - _C*(λa*x_0ij_λa*(aS₁_a + B_a) - λr*x_0ij_λr*(aS₁_r + B_r))
        θ = expm1(τ)
        γc = 10 * (-tanh(10*(0.57 - α)) + 1) * _ζst*θ*exp(_ζst*(-6.7 - 8*_ζst))
        ∂a_2∂ρ_S = 0.5*_C^2 *
                (ρS_∂KHS*(x_0ij_2λa*(aS₁_2a+B_2a)
                        - 2*x_0ij_λaλr*(aS₁_ar+B_ar)
                        + x_0ij_2λr*(aS₁_2r+B_2r)
                    )
                + _KHS*(x_0ij_2λa*(∂aS₁∂ρS_2a + ∂B∂ρS_2a)
                        - 2*x_0ij_λaλr*(∂aS₁∂ρS_ar + ∂B∂ρS_ar)
                        + x_0ij_2λr*(∂aS₁∂ρS_2r + ∂B∂ρS_2r)
                    )
                )

        
        gMCA2 = 3*∂a_2∂ρ_S-_KHS*_C^2 *
                (λr*x_0ij_2λr*(aS₁_2r+B_2r) -
                    (λa+λr)*x_0ij_λaλr*(aS₁_ar+B_ar) +
                    λa*x_0ij_2λa*(aS₁_2a+B_2a)
                )
        g_2_ = (1 + γc)*gMCA2
        log_g_Mie_ = log(g_HSi)+(τ*g_1_/g_HSi+τ^2*g_2_/g_HSi)
        achain -= z[i]*(log_g_Mie_*(mi - 1))
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
            (x_0ij^λa*(@f(aS_1,λa,ζₓ)+@f(B,λa,x_0ij,ζₓ)) - x_0ij^λr*(@f(aS_1,λr,ζₓ)+@f(B,λr,x_0ij,ζₓ)))

            #calculations for a2
            α = _C*(1/(λa-3)-1/(λr-3))
            f1,f2,f3,f4,f5,f6 = @f(f123456,α)
            _χ = f1*_ζst+f2*_ζst5+f3*_ζst8
            a2_ij = π*_KHS*(1+_χ)*ρS*ϵ^2*dij3*_C^2 *
            (x_0ij^(2*λa)*(@f(aS_1,2*λa,ζₓ)+@f(B,2*λa,x_0ij,ζₓ))
            - 2*x_0ij^(λa+λr)*(@f(aS_1,λa+λr,ζₓ)+@f(B,λa+λr,x_0ij,ζₓ))
            + x_0ij^(2*λr)*(@f(aS_1,2λr,ζₓ)+@f(B,2*λr,x_0ij,ζₓ)))

            #calculations for a3
            a3_ij = -ϵ^3*f4*_ζst * exp(_ζst*(f5+f6*_ζst))
            #adding
            a₁ += 2*a1_ij*x_Si*x_Sj
            a₂ += 2*a2_ij*x_Si*x_Sj
            a₃ += 2*a3_ij*x_Si*x_Sj
        end
    end
    a₁ = a₁*m̄/T/∑z
    a₂ = a₂*m̄/(T*T)/∑z
    a₃ = a₃*m̄/(T*T*T)/∑z
    adisp = a₁ + a₂ + a₃
    return adisp + achain/∑z
end

function a_disp(model::SAFTVRMieModel, V, T, z,_data = @f(data))
    _d,ρS,ζi,_ζ_X,_ζst,_,m̄ = _data
    comps = 1:length(z)
    #this is a magic trick. we normally (should) expect length(z) = length(model),
    #but on GC models, @comps != @groups
    #if we pass Xgc instead of z, the equation is exactly the same.
    #we need to add the divide the result by sum(z) later.
    m = model.params.segment.values
    _ϵ = model.params.epsilon.values
    _λr = model.params.lambda_r.values
    _λa = model.params.lambda_a.values
    _σ = model.params.sigma.values
    m̄inv = 1/m̄
    a₁ = zero(V+T+first(z)+one(eltype(model)))
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
    adisp = a₁ + a₂ + a₃
    return adisp
end

function a_chain(model::SAFTVRMieModel, V, T, z,_data = @f(data))
    _d,ρS,ζi,_ζ_X,_ζst,_,m̄ = _data
    l = length(z)
    comps = 1:l
    ∑z = ∑(z)
    m = model.params.segment.values
    _ϵ = model.params.epsilon.values
    _λr = model.params.lambda_r.values
    _λa = model.params.lambda_a.values
    _σ = model.params.sigma.values
    m̄inv = 1/m̄
    a₁ = zero(V+T+first(z)+one(eltype(model)))
    a₂ = a₁
    a₃ = a₁
    achain = a₁
    _ζst5 = _ζst^5
    _ζst8 = _ζst^8
    _KHS,ρS_∂KHS = @f(KHS_f_ρdf,_ζ_X)
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
        γc = 10 * (-tanh(10*(0.57-α))+1) * _ζst*θ*exp(-6.7*_ζst-8*_ζst^2)
        ∂a_2∂ρ_S = 0.5*_C^2 *
            (ρS_∂KHS*(x_0ij^(2*λa)*(aS_1_2a+B_2a)
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
        log_g_Mie_ = log(g_HSi) + (ϵ/T*g_1_/g_HSi+(ϵ/T)^2*g_2_/g_HSi)
        #@show (g_Mie_,i)
        achain +=  z[i]*(log_g_Mie_*(m[i]-1))
    end
    return -achain/∑z
end

########
#=
Optimizations for single component SAFTVRMie
=#

#######

function d(model::SAFTVRMie, V, T, z::SingleComp)
    ϵ = model.params.epsilon.values[1,1]
    σ = model.params.sigma.values[1,1]
    λa = model.params.lambda_a.values[1,1]
    λr = model.params.lambda_r.values[1,1]
    return SA[d_vrmie(T,λa[1],λr[1],σ[1],ϵ[1])]
end
