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
default_references(::Type{CPPCSAFT}) = ["10.1021/ie0003887", "10.1021/ie010954d"]
default_locations(::Type{CPPCSAFT}) = ["SAFT/PCSAFT/CPPCSAFT/", "properties/molarmass.csv", "properties/critical.csv"]
function transform_params(::Type{CPPCSAFT},params,components)
    k = get(params,"k",nothing)
    l = get(params,"l",nothing)
    sigma,epsilon,segment = params["sigma"],params["epsilon"],params["segment"]
    params["sigma"] = sigma_LorentzBerthelot(sigma)
    params["segment"] = sigma_LorentzBerthelot(segment, l)
    params["epsilon"] = epsilon_LorentzBerthelot(epsilon, k)
    return params
end

function a_res(model::CPPCSAFTModel, V, T, z)
    _data = @f(data)
    return @f(a_hc,_data) + @f(a_disp,_data) + @f(a_assoc,_data)
end

function data(model::CPPCSAFTModel, V, T, z)
    m = model.params.segment
    ϵ = model.params.epsilon.values
    σ = model.params.sigma.values
    m̄ = dot(z,diagvalues(m))
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
    (_d,ζ0,ζ1,ζ2,ζ3,m̄,ϵmix,σmix)
end

function a_hs(model::CPPCSAFTModel,V,T,z,_data = @f(data))
    _d,ζ0,ζ1,ζ2,ζ3,m̄,ϵmix,σmix = _data
    θ = CPPCSAFT_theta(T,ϵmix)
    onemζ3 = 1 - ζ3
    RT = Rgas(model)*T
    Ahs = (
        RT * m̄ / ζ0
        * (
            3 * ζ1 * ζ2 / onemζ3
            + ζ2^3 / (ζ3 * onemζ3^2)
            + (ζ2^3 / ζ3^2 - ζ₀) * log1p(-ζ3)
        ) * sqrt(onemζ3 / (1 - ζ3 / θ^3))
    )
    return Ahs
end

function a_hc(model::PCSAFTModel, V, T, z,_data=@f(data))
    dii,ζ0,ζ1,ζ2,ζ3,m̄,_,_ = _data
    m = model.params.segment.values
    Σz = sum(z)
    c1 = 1/(1-ζ3)
    c2 = 3ζ2/(1-ζ3)^2
    c3 = 2ζ2^2/(1-ζ3)^3
    a_hs = @f(a_hs,_data)
    res = zero(a_hs)
    for i ∈ @comps
        dᵢ = dii[i]
        zᵢ = z[i]
        xdᵢ = dᵢ/2
        g_hsᵢᵢ = c1 + c2*xdᵢ + c3*xdᵢ*xdᵢ
        res += zᵢ*zᵢ*(m[i,i]-1)*log(g_hsᵢᵢ)
        for j ∈ 1:(i-1)
            zⱼ = z[j]
            dⱼ = dii[j]
            xdᵢⱼ = dᵢ*dⱼ/(dᵢ+dⱼ)
            g_hsᵢⱼ = c1 + c2*xdᵢⱼ + c3*xdᵢ*xdᵢⱼ
            res += zᵢ*zⱼ*(m[i,j]-1)*log(g_hsᵢⱼ)
        end
    end
    #return  m̄*@f(a_hs) - ∑(z[i]*(m[i]-1)*log(@f(g_hs,i,i)) for i ∈ @comps)/Σz
    return a_hs/Σz - res/(Σz*Σz)
end

function a_disp(model::CPPCSAFTModel, V, T, z,_data=@f(data))
    dii,ζ0,ζ1,ζ2,ζ3,m̄,ϵ,σ = _data
    Σz = sum(z)
    m2ϵσ3 = m̄*m̄*ϵ*σ*σ*σ
    πNAρ = π*N_A*Σz/V
    C₁ = C1(model, V, T, z, ζ3, m̄)
    return -2*πNAρ*@f(I,1,_data)*m2ϵσ3 - m̄*πNAρ*C₁*@f(I,2,_data)*m2ϵσ3/T
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
    ϵᵢᵢ = diagvalues(model.params.epsilon)[1]
    σᵢᵢ = diagvalues(model.params.sigma)[1]
    return SA[σᵢᵢ*CPPCSAFT_theta(T,ϵᵢᵢ)]
end

function CPPCSAFT_theta(T,ϵ)
    Tr = T/ϵ
    (1 + 0.2977 * Tr) / (1 + Tr * (0.33163 + Tr * 0.0010477))
end

function I(model::PCSAFTModel, V, T, z, n, _data=@f(data))
    dii,ζ0,ζ1,ζ2,ζ3,m̄,_,_ = _data
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
