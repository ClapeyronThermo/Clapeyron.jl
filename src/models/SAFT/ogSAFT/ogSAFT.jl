struct ogSAFTParam <: EoSParam
    Mw::SingleParam{Float64}
    segment::SingleParam{Float64}
    sigma::PairParam{Float64}
    epsilon::PairParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
end

abstract type ogSAFTModel <: SAFTModel end
@newmodel ogSAFT ogSAFTModel ogSAFTParam
default_references(::Type{ogSAFT}) = ["10.1021/ie00104a021","10.1016/0378-3812(89)80308-5"]
default_locations(::Type{ogSAFT}) = ["SAFT/ogSAFT","properties/molarmass.csv"]
function transform_params(::Type{ogSAFT},params)
    sigma = params["sigma"]
    sigma.values .*= 1E-10
    return saft_lorentz_berthelot(params)
end

function get_k(model::ogSAFT)   
    return get_k_geomean(model.params.epsilon)
end

function get_l(model::ogSAFT)   
    return get_k_mean(model.params.sigma)
end

"""
    ogSAFTModel <: SAFTModel

    ogSAFT(components;
    idealmodel = BasicIdeal,
    userlocations = String[],
    ideal_userlocations = String[],
    reference_state = nothing,
    verbose = false,
    assoc_options = AssocOptions())

## Input parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Single Parameter (`Float64`) - Segment Diameter [`A°`]
- `epsilon`: Single Parameter (`Float64`) - Reduced dispersion energy  `[K]`
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Paramater (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m^3]`

## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Pair Parameter (`Float64`) - Mixed segment Diameter `[m]`
- `epsilon`: Pair Parameter (`Float64`) - Mixed reduced dispersion energy`[K]`
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume

## Input models
- `idealmodel`: Ideal Model

## Description

(original) Statistical Associating Fluid Theory (og-SAFT) Equation of State

## References
1. Chapman, W. G., Gubbins, K. E., Jackson, G., & Radosz, M. (1989). SAFT: Equation-of-state solution model for associating fluids. Fluid Phase Equilibria, 52, 31–38. [doi:10.1016/0378-3812(89)80308-5](https://doi.org/10.1016/0378-3812(89)80308-5)
2. Chapman, W. G., Gubbins, K. E., Jackson, G., & Radosz, M. (1990). New reference equation of state for associating liquids. Industrial & Engineering Chemistry Research, 29(8), 1709–1721. [doi:10.1021/ie00104a021](https://doi.org/10.1021/ie00104a021)
"""
ogSAFT

export ogSAFT

recombine_impl!(model::ogSAFTModel) = recombine_saft!(model)

function data(model::ogSAFTModel, V, T, z)
    _d = @f(d)
    m̄ = dot(z,model.params.segment.values)
    ζi = @f(ζ0123,_d)
    return _d, m̄, ζi
end

function a_res(model::ogSAFTModel, V, T, z, _data = @f(data))
    return @f(a_seg,_data) + @f(a_chain,_data) + @f(a_assoc,_data)
end

function a_seg(model::ogSAFTModel, V, T, z,_data = @f(data))
    _d, m̄, ζi = _data
    return m̄*(@f(a_hs,_data)+@f(a_disp,_data))/sum(z)
end

function a_chain(model::ogSAFTModel, V, T, z)
    #x = z/∑(z)
    m = model.params.segment.values
    return sum(z[i]*(1-m[i])*log(@f(g_hsij, i, i)) for i ∈ @comps)/sum(z)
end

function d(model::ogSAFTModel, V, T, z)
    ϵ = model.params.epsilon.values
    σ = model.params.sigma.values
    di = zeros(eltype(T+one(eltype(model))),length(model))
    for i in 1:length(model)
        di[i] = σ[i,i]*@f(f_d,ϵ[i,i])
    end
    return di
end

function d(model::ogSAFT, V, T, z::SingleComp)
    ϵ = only(model.params.epsilon.values)
    σ = only(model.params.sigma.values)
    return SA[σ*@f(f_d,ϵ)]
end

function f_d(model::ogSAFTModel, V, T, z,ϵ = model.params.epsilon.values[i,i])
    fm = 0.0010477#+0.025337*(m[i]-1)/m[i]
    f = (1+0.2977T/ϵ)/(1+0.33163T/ϵ+fm*(T/ϵ)^2)
    return f
end

function ζ0123(model::ogSAFTModel, V, T, z,_d=@f(d))
    #N_A*π/6/V * sum(z[i]*m[i]*@f(d,i)^n for i ∈ @comps)
    m = model.params.segment
    _0 = zero(V+T+first(z)+one(eltype(model)))
    ζ0,ζ1,ζ2,ζ3 = _0,_0,_0,_0
    for i ∈ 1:length(z)
        di =_d[i]
        xS = z[i]*m[i]
        ζ0 += xS
        ζ1 += xS*di
        ζ2 += xS*di*di
        ζ3 += xS*di*di*di
    end
    c = π/6*N_A/V
    ζ0,ζ1,ζ2,ζ3 = c*ζ0,c*ζ1,c*ζ2,c*ζ3
    return ζ0,ζ1,ζ2,ζ3
end


# function dx(model::ogSAFTModel, V, T, z)
#     x = z/∑(z)
#     m = model.params.segment.values
#     σ = model.params.sigma.values
#     ϵ = model.params.epsilon.values
#     comps = @comps
#     mx = ∑(x .* m)
#     σx = (∑(x[i]*x[j]*m[i]*m[j]*σ[i,j]^3 for i ∈ comps for j ∈ comps)/mx^2)^(1/3)
#     ϵx = (∑(x[i]*x[j]*m[i]*m[j]*σ[i,j]^3*ϵ[i,j] for i ∈ comps for j ∈ comps)/mx^2)/σx^3

#     fm = 0.0010477#+0.025337*(mx-1)/mx
#     f = (1+0.2977T/ϵx)/(1+0.33163T/ϵx+fm*(T/ϵx)^2)
#     return σx * f
# end

function ζn(model::ogSAFTModel, V, T, z, n)
    #∑z = ∑(z)
    #x = z/∑z
    m = model.params.segment.values
    return N_A*π/6/V * ∑(z[i]*m[i]*@f(d, i)^n for i ∈ @comps)
end

# function η(model::ogSAFTModel, V, T, z)
#     ∑z = ∑(z)
#     x = z/∑z
#     m = model.params.segment.values
#     m̄ = ∑(x .* m)
#     return N_A*∑z*π/6/V*@f(dx)^3*m̄
# end

function g_hsij(model::ogSAFTModel, V, T, z, i, j,_data = @f(data))
    _d, m̄, ζi = _data
    ζ0,ζ1,ζ2,ζ3 = ζi
    di,dj = _d[i],_d[j]
    return 1/(1-ζ3) + di*dj/(di+dj)*3ζ2/(1-ζ3)^2 + (di*dj/(di+dj))^2*2ζ2^2/(1-ζ3)^3
end

# function a_hs(model::ogSAFTModel, V, T, z)
#     ηx = @f(η)
#     return (4ηx-3ηx^2)/(1-ηx)^2
# end

function a_hs(model::ogSAFTModel, V, T, z,_data = @f(data))
    _d, m̄, ζi = _data
    ζ0,ζ1,ζ2,ζ3 = ζi
    return bmcs_hs(ζ0,ζ1,ζ2,ζ3)
end

function a_disp(model::ogSAFTModel, V, T, z,_data = @f(data))
    _d, m̄, ζi = _data
    ζ0,ζ1,ζ2,ζ3 = ζi

    m = model.params.segment.values
    σ = model.params.sigma.values
    ϵ = model.params.epsilon.values
    mσϵ3 = zero(first(z)+one(eltype(model)))
    mσ3 = zero(first(z)+one(eltype(model)))
    for i in @comps
        zi,ϵi,σi,mi = z[i],ϵ[i,i],σ[i,i],m[i]
        mσ3ii = zi*zi*mi*σi*σi*σi
        mσ3 += mσ3ii
        mσϵ3 += mσ3ii*ϵi
        for j in 1:(i-1)
            zj,ϵij,σij,mj = z[j],ϵ[i,j],σ[i,j],m[j]
            mσ3ij = zi*zj*mi*mj*σij*σij*σij
            mσ3 += 2*mσ3ij
            mσϵ3 += 2*mσ3ij*ϵij
        end
    end

    ϵx = mσϵ3/mσ3
    ρR = (6/sqrt(2)/π)*ζ3
    TR = T/ϵx
    a_seg1 = ρR*evalpoly(ρR,(-8.5959,-4.5424,-2.1268,10.285))
    a_seg2 = ρR*evalpoly(ρR,(-1.9075,9.9724,-22.216,+15.904))
    return 1/TR*(a_seg1+a_seg2/TR)
end

## This is an attempt to make Twu et al.'s segment term; does not work yet
# function a_seg(model::ogSAFTModel, V, T, z)
#     Bo = [1.31024,-3.80636,-2.37238,-0.798872,0.198761,1.47014,-0.786367,2.19465,5.75429,6.7822,-9.94904,-15.6162,86.643,18.527,9.04755,8.68282]
#     Ba = [3.79621,-6.14518,-1.84061,-2.77584,-0.420751,-5.66128,19.2144,-3.33443,33.0305,-5.90766,9.55619,-197.883,-61.2535,77.1802,-6.57983,0.0]
#     ω = 0.011
#     A = []
#     for i ∈ 1:16
#         append!(A,Bo[i]+ω*Ba[i])
#     end
#     m = model.params.segment
#     σ = model.params.sigma
#     ϵ = model.params.epsilon
#     x = z/∑(z[i] for i ∈ @comps)
#     mx = ∑(x[i]*m[i] for i ∈ @comps)
#     σx = (∑(x[i]*x[j]*m[i]*m[j]*σ[i,j]^3 for i ∈ @comps for j ∈ @comps)/mx^2)^(1/3)
#     ϵx = (∑(x[i]*x[j]*m[i]*m[j]*σ[i,j]^3*ϵ[i,j] for i ∈ @comps for j ∈ @comps)/mx^2)/σx^3
#
#     ρ  = ∑(z)*N_A/V
#     # ρR = ρ*mx*σx^3
#     ηx = η(model,V, T, z)
#     ρR = (6/π)*ηx
#     TR = T/ϵx
#
#     u_res = (A[2]/TR+2A[3]/TR^2+3A[4]/TR^3+5A[5]/TR^5)*ρR+1/2*A[7]/TR*ρR^2+
#             1/(2*A[16])*(3A[9]/TR^3+4A[10]/TR^4+5A[11]/TR^5)*(1-exp(-A[16]*ρR^2))+
#             1/(2*A[16]^2)*(3A[12]/TR^3+4A[13]/TR^4+5A[14]/TR^5)*(1-(1+A[16]*ρR^2)*exp(-A[16]*ρR^2))+
#             1/5*A[15]/TR*ρR^5
#     s_res = -log(ρ*R̄*T)-(A[1]-A[3]/TR^2-2A[4]/TR^3-4A[5]/TR^5)*ρR-1/2*A[6]*ρR^2-1/3*A[8]*ρR^3+
#             1/(2*A[16])*(2A[9]/TR^3+3A[10]/TR^4+4A[11]/TR^5)*(1-exp(-A[16]*ρR^2))+
#             1/(2*A[16]^2)*(2A[12]/TR^3+3A[13]/TR^4+4A[14]/TR^5)*(1-(1+A[16]*ρR^2)*exp(-A[16]*ρR^2))
#     a_res = u_res-s_res
#     return mx*(a_res)
# end

function Δ(model::ogSAFTModel, V, T, z, i, j, a, b,_data = @f(data))
    _d, m̄, ζi = _data
    ζ0,ζ1,ζ2,ζ3 = ζi
    κ = model.params.bondvol.values
    kijab =κ[i,j][a,b]
    ϵ_assoc = model.params.epsilon_assoc.values
    g = @f(g_hsij,i,j)
    return (_d[i]+_d[j])^3/2^3*g*(expm1(ϵ_assoc[i,j][a,b]/T))*kijab
end

