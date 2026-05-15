struct PCSAFTParam{T} <: ParametricEoSParam{T}
    Mw::SingleParam{T}
    segment::SingleParam{T}
    sigma::PairParam{T}
    epsilon::PairParam{T}
    epsilon_assoc::AssocParam{T}
    bondvol::AssocParam{T}
end

function PCSAFTParam(Mw,segment,sigma,epsilon,epsilon_assoc,bondvol)
    return build_parametric_param(PCSAFTParam,Mw,segment,sigma,epsilon,epsilon_assoc,bondvol)
end

abstract type PCSAFTModel <: SAFTModel end
@newmodel PCSAFT PCSAFTModel PCSAFTParam{T}
default_references(::Type{PCSAFT}) = ["10.1021/ie0003887", "10.1021/ie010954d"]
default_locations(::Type{PCSAFT}) = ["SAFT/PCSAFT","properties/molarmass.csv"]

function transform_params(::Type{PCSAFT},params)
    sigma = params["sigma"]
    sigma.values .*= 1E-10
    return saft_lorentz_berthelot(params)
end

function x0_volume_liquid(model::PCSAFTModel,p,T,z)
    lb_v = lb_volume(model,T,z)
    Ts = T_scale(model,z)
    if T > 0.9Ts
        return 1.25*lb_v
    else
        return x0_volume_liquid_lowT(model,p,T,z)
    end
end

function get_k(model::PCSAFTModel)
    has_groups(model) && return nothing
    return get_k_geomean(model.params.epsilon)
end

function get_l(model::PCSAFTModel)
    has_groups(model) && return nothing
    return get_k_mean(model.params.sigma)
end

function set_k!(model::PCSAFTModel,k)
    has_groups(model) && return nothing
    return recombine_saft!(model,k,nothing)
end

function set_l!(model::PCSAFTModel,l)
    has_groups(model) && return nothing
    return recombine_saft!(model,nothing,l)
end


"""
    PCSAFTModel <: SAFTModel

    PCSAFT(components;
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
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Parameter (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m³]`

## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Pair Parameter (`Float64`) - Mixed segment Diameter `[m]`
- `epsilon`: Pair Parameter (`Float64`) - Mixed reduced dispersion energy `[K]`
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m³]`

## Input models
- `idealmodel`: Ideal Model

## Description
Perturbed-Chain SAFT (PC-SAFT)

## References
1. Gross, J., & Sadowski, G. (2001). Perturbed-chain SAFT: An equation of state based on a perturbation theory for chain molecules. Industrial & Engineering Chemistry Research, 40(4), 1244–1260. [doi:10.1021/ie0003887](https://doi.org/10.1021/ie0003887)
2. Gross, J., & Sadowski, G. (2002). Application of the perturbed-chain SAFT equation of state to associating systems. Industrial & Engineering Chemistry Research, 41(22), 5510–5515. [doi:10.1021/ie010954d](https://doi.org/10.1021/ie010954d)
"""
PCSAFT

export PCSAFT

recombine_impl!(model::PCSAFTModel) = recombine_saft!(model)

function a_res(model::PCSAFTModel, V, T, z, _data = @f(data))
    return @f(a_hc,_data) + @f(a_disp,_data) + @f(a_assoc,_data)
end

function data(model::PCSAFTModel,V,T,z)
    _d = @f(d)
    ζ0,ζ1,ζ2,ζ3 = @f(ζ0123,_d)
    m = model.params.segment.values
    m̄ = dot(z, m)/sum(z)
    return (_d,ζ0,ζ1,ζ2,ζ3,m̄)
end

#unpacks packing_fraction from data
function packing_fraction(model::PCSAFTModel,_data)
    _,_,_,_,η,_ = _data
    return η
end

function a_hc(model::PCSAFTModel, V, T, z,_data=@f(data))
    dii,ζ0,ζ1,ζ2,ζ3,m̄ = _data
    m = model.params.segment.values
    Σz = sum(z)
    c1 = 1/(1-ζ3)
    c2 = 3ζ2/(1-ζ3)^2
    c3 = 2ζ2^2/(1-ζ3)^3
    if !iszero(ζ3)
        a_hs = bmcs_hs(ζ0,ζ1,ζ2,ζ3)
    else
        a_hs = @f(bmcs_hs_zero_v,dii)
    end
    res = zero(a_hs)
    for i ∈ @comps
        dᵢ = dii[i]
        di,dj = dᵢ,dᵢ
        g_hsᵢᵢ = c1 + di*dj/(di+dj)*c2 + (di*dj/(di+dj))^2*c3
        res += z[i]*(m[i]-1)*log(g_hsᵢᵢ)
    end
    #return  m̄*@f(a_hs) - ∑(z[i]*(m[i]-1)*log(@f(g_hs,i,i)) for i ∈ @comps)/Σz
    return m̄*a_hs - res/Σz
end

function a_disp(model::PCSAFTModel, V, T, z,_data=@f(data))
    di,ζ0,ζ1,ζ2,ζ3,m̄ = _data
    Σz = sum(z)
    m2ϵσ3₁,m2ϵσ3₂ = @f(m2ϵσ3)
    πNAρ = π*N_A*Σz/V
    return -2*πNAρ*@f(I,1,_data)*m2ϵσ3₁ - m̄*πNAρ*@f(C1,_data)*@f(I,2,_data)*m2ϵσ3₂
end

d(model::PCSAFTModel, V, T, z) = ck_diameter(model, T, z)

#defined in SAFT/equations.jl
#=
function ζ(model::PCSAFTModel, V, T, z, n, _d = @f(d))
    m = model.params.segment.values
    res = zero(V+T+first(z)+one(eltype(model)))
    for i ∈ @comps
        dᵢ = _d[i]
        res += z[i]*m[i]*dᵢ^n
    end
    res *= N_A*π/6/V
    return res
end =#

function g_hs(model::PCSAFTModel, V, T, z, i, j, _data=@f(data))
    _d,ζ0,ζ1,ζ2,ζ3,_ = _data
    return g_hs_ij(_d,ζ2,ζ3,i,j)
end

function a_hs(model::PCSAFTModel, V, T, z, _data=@f(data))
    _,ζ0,ζ1,ζ2,ζ3,_ = _data
    return bmcs_hs(ζ0,ζ1,ζ2,ζ3)
end

function C1(model::PCSAFTModel, V, T, z, _data=@f(data))
    _,_,_,_,η,m̄ = _data
    return C1(model, V, T, z, η, m̄)
end

function C1(model::PCSAFTModel, V, T, z, η, m̄)
    return (1 + m̄*(8η-2η^2)/(1-η)^4 + (1-m̄)*evalpoly(η,(0,20,-27,12,-2))/((1-η)*(2-η))^2)^-1
end

function m2ϵσ3(model::PCSAFTModel, V, T, z)
    m = model.params.segment.values
    σ = model.params.sigma.values
    ϵ = model.params.epsilon.values
    m2ϵσ3₂ = zero(T+first(z)+one(eltype(model)))
    m2ϵσ3₁ = m2ϵσ3₂
    @inbounds for i ∈ @comps
        for j ∈ @comps
            constant = z[i]*z[j]*m[i]*m[j] * σ[i,j]^3
            exp1 = (ϵ[i,j]/T)
            exp2 = exp1*exp1
            m2ϵσ3₁ += constant*exp1
            m2ϵσ3₂ += constant*exp2
        end
    end
    Σz = sum(z)
    k = (1/Σz)^2
    return k*m2ϵσ3₁,k*m2ϵσ3₂
    #return ∑(z[i]*z[j]*m[i]*m[j] * (ϵ[i,j]*(1)/T)^n * σ[i,j]^3 for i ∈ @comps, j ∈ @comps)/(sum(z)^2)
end

function I(model::PCSAFTModel, V, T, z, n, _data=@f(data))
    _,_,_,_,η,m̄ = _data
    if n == 1
        corr = PCSAFTconsts.corr1
    elseif n == 2
        corr = PCSAFTconsts.corr2
    end
    res = zero(η)
    m̄1 = (m̄-1)/m̄
    m̄2 = (m̄-1)/m̄*(m̄-2)/m̄
    ηi = one(η)
    @inbounds for i ∈ 1:length(corr)
        corr1,corr2,corr3 = corr[i]
        ki = corr1 + m̄1*corr2 + m̄2*corr3
        res += ki*ηi
        ηi *= η
    end
    return res
end

function Δ(model::PCSAFTModel, V, T, z, i, j, a, b,_data=@f(data))
    _0 = zero(V+T+first(z)+one(eltype(model)))
    ϵ_assoc = model.params.epsilon_assoc.values
    κ = model.params.bondvol.values
    κijab = κ[i,j][a,b]
    iszero(κijab) && return _0
    σ = model.params.sigma.values
    gij = @f(g_hs,i,j,_data)
    res = gij*σ[i,j]^3*(expm1(ϵ_assoc[i,j][a,b]/T))*κijab
    return res
end

const PCSAFTconsts = (
    corr1 =
    SA[(0.9105631445,-0.3084016918, -0.0906148351),
    (0.6361281449, 0.1860531159, 0.4527842806),
    (2.6861347891, -2.5030047259, 0.5962700728),
    (-26.547362491, 21.419793629, -1.7241829131),
    (97.759208784, -65.255885330, -4.1302112531),
    (-159.59154087, 83.318680481, 13.776631870),
    (91.297774084, -33.746922930, -8.6728470368)],

    corr2 =
    SA[(0.7240946941, -0.5755498075, 0.0976883116),
    (2.2382791861, 0.6995095521, -0.2557574982),
    (-4.0025849485, 3.8925673390, -9.1558561530),
    (-21.003576815, -17.215471648, 20.642075974),
    (26.855641363, 192.67226447, -38.804430052),
    (206.55133841, -161.82646165, 93.626774077),
    (-355.60235612, -165.20769346, -29.666905585)]
)

#=
Specific PCSAFT optimizations
This code is not generic, in the sense that is only used by PCSAFT and not any model <:PCSAFTModel
but, because it is one of the more commonly used EoS,
It can have some specific optimizations to make it faster.
=#

#Optimized Δ for PCSAFT

function  Δ(model::PCSAFT, V, T, z,_data=@f(data))
    ϵ_assoc = model.params.epsilon_assoc.values
    κ = model.params.bondvol.values
    σ = model.params.sigma.values
    Δout = assoc_similar(κ,typeof(V+T+first(z)+one(eltype(model))))
    Δout.values .= false  #fill with zeros, maybe it is not necessary?
    for (idx,(i,j),(a,b)) in indices(Δout)
        κijab = κ[idx]
        if κijab != 0
            gij = @f(g_hs,i,j,_data)
            Δout[idx] = gij*σ[i,j]^3*(expm1(ϵ_assoc[i,j][a,b]/T))*κijab
        end
    end
    return Δout
end
