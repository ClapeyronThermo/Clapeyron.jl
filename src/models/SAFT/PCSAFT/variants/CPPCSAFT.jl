struct CPPCSAFTParam <: EoSParam
    Mw::SingleParam{Float64}
    segment::PairParam{Float64}
    sigma::PairParam{Float64}
    epsilon::PairParam{Float64}
    delta::SingleParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
end

abstract type CPPCSAFTModel <: PCSAFTModel end
@newmodel CPPCSAFT CPPCSAFTModel CPPCSAFTParam
export CPPCSAFT

"""
    CPPCSAFTModel <: PCSAFTModel

    CPPCSAFT(components; 
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
- `delta`: Single Parameter (`Float64`) - Critical volume displacement (no units)
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Paramater (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m^3]`
## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `segment`: Pair Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Pair Parameter (`Float64`) - Mixed segment Diameter `[m]`
- `epsilon`: Pair Parameter (`Float64`) - Mixed reduced dispersion energy`[K]`
- `delta`: Single Parameter (`Float64`) - Critical volume displacement (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume
## Input models
- `idealmodel`: Ideal Model
## Description
Critical Point-Based PC-SAFT (CP-PC-SAFT)
## References
1. Polishuk, I. (2014). Standardized critical point-based numerical solution of statistical association fluid theory parameters: The perturbed chain-statistical association fluid theory equation of state revisited. Industrial & Engineering Chemistry Research, 53(36), 14127–14141. [doi:10.1021/ie502633e](https://doi.org/10.1021/ie502633e)
"""
CPPCSAFT
default_references(::Type{CPPCSAFT}) = ["10.1021/ie502633e"]
default_locations(::Type{CPPCSAFT}) = ["SAFT/PCSAFT/CPPCSAFT/", "properties/molarmass.csv"]
function transform_params(::Type{CPPCSAFT},params,components)
    k = get(params,"k",nothing)
    l = get(params,"l",nothing)
    sigma,epsilon,segment = params["sigma"],params["epsilon"],params["segment"]
    sigma.values .*= 1E-10
    params["sigma"] = sigma_LorentzBerthelot(sigma)
    params["segment"] = sigma_LorentzBerthelot(segment, l)
    params["epsilon"] = epsilon_LorentzBerthelot(epsilon, k)
    return params
end

function get_l(model::CPPCSAFTModel)   
    return get_k_mean(model.params.segment)
end

function recombine_impl!(model::CPPCSAFTModel)
    sigma = model.params.sigma
    epsilon = model.params.epsilon
    segment = model.params.segment
    sigma = sigma_LorentzBerthelot!(sigma)
    segment = sigma_LorentzBerthelot!(segment)
    epsilon = epsilon_LorentzBerthelot!(epsilon)
    return model
end

function a_res(model::CPPCSAFTModel, V, T, z)
    _data = @f(data)
    return @f(a_hc,_data) + @f(a_disp,_data) + @f(a_assoc,_data)
end

function data(model::CPPCSAFTModel, V, T, z)
    m = model.params.segment.values
    ϵ = model.params.epsilon.values
    σ = model.params.sigma.values
    m̄ = sum(z[i]*m[i,i] for i in @comps)
    ϵmix = zero(eltype(model))
    σmix = zero(eltype(model))
    for i in @comps
        ϵi,σi,zi,mi = ϵ[i,i],σ[i,i],z[i],m[i,i]
        σmix_ii = zi*zi*mi*mi*(σi*σi*σi)
        ϵmix += σmix_ii*ϵi
        σmix += σmix_ii
        for j in 1:(i-1)
            ϵij,σij,zj,mj = ϵ[i,j],σ[i,j],z[j],m[j,j]
            σmix_ij = zi*zj*mi*mj*(σij*σij*σij)
            ϵmix += 2*zi*zj*mi*mj*(σij*σij*σij)*ϵij
            σmix += 2*σmix_ij
        end
    end
    σmix = cbrt(σmix/(m̄*m̄)) #no units
    ϵmix = ϵmix/(σmix*σmix*σmix*m̄*m̄)
    m̄ = m̄/sum(z)
    _d = @f(d)
    ζ0,ζ1,ζ2,ζ3 = @f(ζ0123,_d,diagvalues(m))
    (_d,ζ0,ζ1,ζ2,ζ3,m̄,ϵmix,σmix)
end

function a_hs(model::CPPCSAFTModel,V,T,z,_data = @f(data))
    _d,ζ0,ζ1,ζ2,ζ3,m̄,ϵmix,σmix = _data
    θ = CPPCSAFT_theta(T,ϵmix)
    _a_hs = bmcs_hs(ζ0,ζ1,ζ2,ζ3)
    return _a_hs*sqrt((1 - ζ3) / (1 - ζ3 / θ^3))
end

function a_hc(model::CPPCSAFTModel, V, T, z,_data=@f(data))
    dii,ζ0,ζ1,ζ2,ζ3,m̄,_,_ = _data
    m = model.params.segment.values
    Σz = sum(z)
    c1 = 1/(1-ζ3)
    c2 = 3ζ2/(1-ζ3)^2
    c3 = 2ζ2^2/(1-ζ3)^3
    _a_hs = @f(a_hs,_data)
    _a_disp = zero(_a_hs)
    for i ∈ @comps
        dᵢ = dii[i]
        zᵢ = z[i]
        xdᵢ = dᵢ/2
        g_hsᵢᵢ = c1 + c2*xdᵢ + c3*xdᵢ*xdᵢ
        _a_disp += zᵢ*zᵢ*(m[i,i]-1)*log(g_hsᵢᵢ)
        for j ∈ 1:(i-1)
            zⱼ = z[j]
            dⱼ = dii[j]
            xdᵢⱼ = dᵢ*dⱼ/(dᵢ+dⱼ)
            g_hsᵢⱼ = c1 + c2*xdᵢⱼ + c3*xdᵢ*xdᵢⱼ
            _a_disp += zᵢ*zⱼ*(m[i,j]-1)*log(g_hsᵢⱼ)
        end
    end
    #return  m̄*@f(a_hs) - ∑(z[i]*(m[i]-1)*log(@f(g_hs,i,i)) for i ∈ @comps)/Σz
    return m̄*_a_hs - _a_disp/(Σz*Σz)
end

function a_disp(model::CPPCSAFTModel, V, T, z,_data=@f(data))
    dii,ζ0,ζ1,ζ2,ζ3,m̄,ϵ,σ = _data
    Σz = sum(z)
    m2ϵσ3 = m̄*m̄*ϵ*σ*σ*σ
    πNAρ = π*N_A*Σz/V
    C₁ = C1(model, V, T, z, ζ3, m̄)
    return -πNAρ*m2ϵσ3*(2*@f(I,1,_data)/T + C₁*m̄*ϵ*@f(I,2,_data)/(T*T))
end

function d(model::CPPCSAFTModel, V, T, z)
    ϵᵢᵢ = diagvalues(model.params.epsilon)
    σᵢᵢ = diagvalues(model.params.sigma)
    di = zeros(eltype(T+one(eltype(model))),length(model))
    for i in eachindex(di)
        di[i] = σᵢᵢ[i]*CPPCSAFT_theta(T,ϵᵢᵢ[i])
    end
    return di
end

function d(model::CPPCSAFTModel, V, T, z::SingleComp)
    ϵ = only(model.params.epsilon.values)
    σ = only(model.params.sigma.values)
    return SA[σ*CPPCSAFT_theta(T,ϵ)]
end

function CPPCSAFT_theta(T,ϵ)
    Tr = T/ϵ
    (1 + 0.2977 * Tr) / (1 + Tr * (0.33163 + Tr * 0.0010477))
end

function I(model::CPPCSAFTModel, V, T, z, n, _data=@f(data))
    dii,ζ0,ζ1,ζ2,η,m̄,_,_ = _data
    if n == 1
        corr = CPPCSAFTconsts.corr1
    elseif n == 2
        corr = CPPCSAFTconsts.corr2
    end
    res = zero(η)
    m̄1 = (m̄-1)/m̄
    m̄2 = (m̄-1)/m̄*(m̄-2)/m̄
    @inbounds for i ∈ 1:length(corr)
        ii = i-1
        corr1,corr2,corr3 = corr[i]
        ki = corr1 + m̄1*corr2 + m̄2*corr3
        res +=ki*η^ii
    end
    return res
end

function lb_volume(model::CPPCSAFTModel, z = SA[1.0])
    m = model.params.segment.values
    σ = model.params.sigma.values
    val = π/6*N_A*sum(z[i]*m[i,i]*σ[i,i]^3 for i in 1:length(z))
    return val
end

const CPPCSAFTconsts = (
    corr1 =
    SA[(0.880823927666, -0.349731891574, -0.041574194083),
    (1.26235042398, 1.06133747189, -0.828880456022),
    (-2.88916037036, -9.92662697237, 10.6610090572),
    (-0.791682734039, 55.1147516007, -42.2676046130),
    (31.4414035626, -158.619888888, 93.3498157944),
    (-67.7739765931, 237.469601780, -119.982855050),
    (37.6471023573, -146.917589624, 69.3982688833)],

    corr2 =
    SA[(0.7240946941, -0.5755498075, 0.0976883116),
    (2.2382791861, 0.6995095521, -0.2557574982),
    (-4.0025849485, 3.8925673390, -9.1558561530),
    (-21.003576815, -17.215471648, 20.642075974),
    (26.855641363, 192.67226447, -38.804430052),
    (206.55133841, -161.82646165, 93.626774077),
    (-355.60235612, -165.20769346, -29.666905585)]
)

function Δ(model::CPPCSAFTModel, V, T, z, i, j, a, b,_data=@f(data))
    _0 = zero(V+T+first(z)+one(eltype(model)))
    _d = first(_data)
    ϵ_assoc = model.params.epsilon_assoc.values
    κ = model.params.bondvol.values
    κijab = κ[i,j][a,b] 
    σ = model.params.sigma.values
    ϵ = model.params.epsilon.values
    iszero(κijab) && return _0
    if i == j #use precalculated values
        dij = _d[i]
    else
        dij = σ[i,j]*CPPCSAFT_theta(T,ϵ[i,j])
    end
    gij = @f(g_hs,i,j,_data)
    res = gij*dij^3*(expm1(ϵ_assoc[i,j][a,b]/T))*κijab
    return res
end

function  Δ(model::CPPCSAFT, V, T, z,_data=@f(data))
    ϵ_assoc = model.params.epsilon_assoc.values
    κ = model.params.bondvol.values
    σ = model.params.sigma.values
    ϵ = model.params.epsilon.values
    _d = first(_data)
    Δout = assoc_similar(κ,typeof(V+T+first(z)+one(eltype(model))))
    Δout.values .= false  #fill with zeros, maybe it is not necessary?
    for (idx,(i,j),(a,b)) in indices(Δout)
        gij = @f(g_hs,i,j,_data)
        if i == j #use precalculated values
            dij = _d[i]
        else
            dij = σ[i,j]*CPPCSAFT_theta(T,ϵ[i,j])
        end
        Δout[idx] = gij*dij^3*(expm1(ϵ_assoc[i,j][a,b]/T))*κ[i,j][a,b]
    end
    return Δout
end
