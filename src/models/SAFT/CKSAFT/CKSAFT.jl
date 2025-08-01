struct CKSAFTParam <: EoSParam
    Mw::SingleParam{Float64}
    segment::SingleParam{Float64}
    sigma::PairParam{Float64}
    epsilon::PairParam{Float64}
    c::SingleParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
end

abstract type CKSAFTModel <: SAFTModel end
@newmodel CKSAFT CKSAFTModel CKSAFTParam
default_references(::Type{CKSAFT}) = ["10.1021/IE00107A014", "10.1021/ie00056a050"]
default_locations(::Type{CKSAFT}) = ["SAFT/CKSAFT","properties/molarmass.csv"]
function transform_params(::Type{CKSAFT},params)
    k = get(params,"k",nothing)
    l = get(params,"l",nothing)
    sigma = params["vol"]
    sigma.values .*= 6*0.74048/N_A/1e6/π
    sigma.values .^= 1/3
    epsilon = params["epsilon"]
    params["sigma"] = sigma_LorentzBerthelot(sigma, l)
    params["epsilon"] = epsilon_LorentzBerthelot(epsilon, k)
    return params
end

function get_k(model::CKSAFT)   
    return get_k_geomean(model.params.epsilon)
end

function get_l(model::CKSAFT)   
    return get_k_mean(model.params.sigma)
end

export CKSAFT
"""
    CKSAFTModel <: SAFTModel

    CKSAFT(components; 
    idealmodel = BasicIdeal,
    userlocations = String[],
    ideal_userlocations = String[],
    reference_state = nothing,
    verbose = false,
    assoc_options = AssocOptions())

## Input parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `vol`: Single Parameter (`Float64`) - Segment Volume `[dm³]`
- `epsilon`: Single Parameter (`Float64`) - Reduced dispersion energy `[K]`
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Paramater (no units)
- `c`: Single Parameter (`Float64`) - Dispersion T-dependent parameter (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m³]`

## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Pair Parameter (`Float64`) - Mixed segment Diameter `[m]`
- `epsilon`: Pair Parameter (`Float64`) - Mixed reduced dispersion energy`[K]`
- `c`: Single Parameter (`Float64`) - Dispersion T-dependent parameter (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m³]`

## Input models
- `idealmodel`: Ideal Model

## Description

Chen and Kreglewski SAFT (CK-SAFT)

## References
1. Huang, S. H., & Radosz, M. (1990). Equation of state for small, large, polydisperse, and associating molecules. Industrial & Engineering Chemistry Research, 29(11), 2284–2294. [doi:10.1021/ie00107a014](https://doi.org/10.1021/ie00107a014)
2. Huang, S. H., & Radosz, M. (1991). Equation of state for small, large, polydisperse, and associating molecules: extension to fluid mixtures. Industrial & Engineering Chemistry Research, 30(8), 1994–2005. [doi:10.1021/ie00056a050](https://doi.org/10.1021/ie00056a050)
"""
CKSAFT

recombine_impl!(model::CKSAFTModel) = recombine_saft!(model)

function data(model::CKSAFTModel, V, T, z)
    _d = @f(d)
    Σz = sum(z)
    m = model.params.segment.values
    m̄ = zero(Σz + one(eltype(model)))
    for i in @comps
        mi,zi = m[i],z[i]
        m̄ += mi*zi*zi
        for j in 1:(i-1)
            mj,zj = m[j],z[j]
            mij = 0.5*(mi + mj)
            m̄ += 2*mij*zi*zj
        end
    end
    m̄ = m̄/Σz/Σz
    ζi = @f(ζ0123,_d)
    return _d, m̄, ζi, Σz
end

function a_res(model::CKSAFTModel, V, T, z, _data = @f(data))
    return @f(a_seg,_data) + @f(a_chain,_data) + @f(a_assoc,_data)
end

d(model::CKSAFTModel, V, T, z) = ck_diameter(model, T, z)

function a_seg(model::CKSAFTModel, V, T, z,_data = @f(data))
    _d, m̄, ζi, Σz = _data
    return m̄*(@f(a_hs,_data) + @f(a_disp,_data))
end

function a_hs(model::CKSAFTModel, V, T, z,_data = @f(data))
    _d, m̄, ζi, Σz = _data
    ζ0,ζ1,ζ2,ζ3 = ζi
    if !iszero(ζ3)
        _a_hs = bmcs_hs(ζ0,ζ1,ζ2,ζ3)
    else
        _a_hs = @f(bmcs_hs_zero_v,_d)
    end
    return _a_hs
end

function a_disp(model::CKSAFTModel, V, T, z,_data = @f(data))
    _d, m̄, ζi, Σz = _data
    ζ0,ζ1,ζ2,ζ3 = ζi
    ϵ̄ = @f(ū,_data)
    η = ζ3
    τ = 0.74048
    D1 = CKSAFT_consts.D1
    D2 = CKSAFT_consts.D2
    D3 = CKSAFT_consts.D3
    D4 = CKSAFT_consts.D4
    ϵT = (ϵ̄/T)
    ητ = η/τ
    A1 = evalpoly(ητ,CKSAFT_consts.D1)
    A0 = zero(A1)
    A2 = evalpoly(ητ,CKSAFT_consts.D2)
    A3 = evalpoly(ητ,CKSAFT_consts.D3)
    A4 = evalpoly(ητ,CKSAFT_consts.D4)
    return evalpoly(ϵT,(A0,A1,A2,A3,A4))
end

function u(model::CKSAFTModel, V, T, z, i, j)
    ϵ0 = model.params.epsilon.values[i,j]
    c = ck_c(model)
    ci,cj,Tinv = c[i],c[j],1/T
    cTi,cTj = 1 + ci*Tinv,1 + cj*Tinv
    if i == j
        return ϵ0*cTi
    else
        return ϵ0*sqrt(cTi*cTj)
    end
end

ck_c(model::CKSAFTModel) = model.params.c.values

function ū(model::CKSAFTModel, V, T, z,_data = @f(data))
    _d, m̄, ζi, Σz = _data
    ζ0,ζ1,ζ2,ζ3 = ζi
    Σz = ∑(z)
    m = model.params.segment.values
    ϵ = model.params.epsilon.values
    c = ck_c(model)
    T⁻¹ = 1/T
    num = zero(V+T+first(z)+one(eltype(model)))
    denom = zero(T+first(z)+one(eltype(model)))
    for i in @comps
        ci,ϵii,mi,zi,di = c[i],ϵ[i,i],m[i],z[i],_d[i]
        cTi = 1 + ci*T⁻¹
        uii = ϵii*cTi
        di3 = di*di*di
        num += zi*zi*mi*mi*uii*di3*di3
        denom += zi*mi*di3
        for j in 1:(i-1)
            cj,ϵij,mj,zj,dj = c[j],ϵ[i,j],m[j],z[j],_d[j]
            cTj = 1 + cj*T⁻¹
            cTij = sqrt(cTi*cTj)
            uij = ϵij*cTij
            dj3 = dj*dj*dj
            num += 2*zi*zj*mi*mj*uij*di3*dj3
        end
    end
    return num/denom/denom
end

function a_chain(model::CKSAFTModel, V, T, z, _data = @f(data))
    _d, m̄, ζi, Σz = _data
    ζ0,ζ1,ζ2,ζ3 = ζi
    m = model.params.segment.values
    return ∑(z[i]*(1-m[i])*log(@f(g_hs,i,i,_data)) for i ∈ @comps)/Σz
end

function g_hs(model::CKSAFTModel, V, T, z, i, j,_data = @f(data))
    _d, m̄, ζi, Σz = _data
    ζ0,ζ1,ζ2,ζ3 = ζi
    return g_hs_ij(_d,ζ2,ζ3,i,j)
end

function Δ(model::CKSAFTModel, V, T, z, i, j, a, b,_data = @f(data))
    ϵ_associjab = model.params.epsilon_assoc.values[i,j][a,b]
    κijab = model.params.bondvol.values[i,j][a,b]
    σij = model.params.sigma.values[i,j]
    gij = @f(g_hs,i,j,_data)
    return gij*σij^3*expm1(ϵ_associjab/T)*κijab
end

const CKSAFT_consts =(
    D1 = [0.0,-8.8043,4.164627,-48.203555,140.4362,-195.23339,113.515],
    D2 = [0.0,2.9396,-6.0865383,40.137956,-76.230797,-133.70055,860.25349,-1535.3224,1221.4261,-409.10539],
    D3 = [0.0,-2.8225,4.7600148,11.257177,-66.382743,69.248785],
    D4 = [0.0,0.34,-3.1875014,12.231796,-12.110681],
    D =
    [0.9105631445 -0.3084016918 -0.0906148351;
    0.6361281449 0.1860531159 0.4527842806;
    2.6861347891 -2.5030047259 0.5962700728;
    -26.547362491 21.419793629 -1.7241829131;
    97.759208784 -65.255885330 -4.1302112531;
    -159.59154087 83.318680481 13.776631870;
    91.297774084 -33.746922930 -8.6728470368],
)
