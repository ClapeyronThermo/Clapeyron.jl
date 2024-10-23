abstract type GrossVrabecQuatrupolarModel <: GrossVrabecModel end

struct GrossVrabecQuatrupolarParam{T} <: ParametricEoSParam{T}
    sigma::PairParam{T}
    epsilon::PairParam{T}
    segment::SingleParam{T}
    dipole::SingleParam{T}
    dipole2::SingleParam{T}
    quadrupole::SingleParam{T}
    quadrupole2::SingleParam{T}
end

@newmodelsimple GrossVrabecQuatrupolar GrossVrabecQuatrupolarModel GrossVrabecQuatrupolarParam

function GrossVrabecQuatrupolar(model::EoSModel)
    components = model.components
    sigma = model.params.sigma
    epsilon = model.params.epsilon
    segment = model.params.segment
    dipole = model.params.dipole
    quatrupole = model.params.quatrupole
    
    paramtype = typeof(model.params)
    if hasfield(paramtype,:dipole2)
        dipole2 = model.params.dipole2
    else
        dipole2 = μ.^2 ./ m ./ k_B*1e-36*(1e-10*1e-3)
    end

    if hasfield(paramtype,:quatrupole2)
        quatrupole2 = model.params.quatrupole2
    else
        quatrupole2 = Q.^2 ./ m ./ k_B*1e-56*(1e-10*1e-3))
    end

    pkg_param = GrossVrabecQuatrupolarParam(sigma,epsilon,segment,dipole,dipole2)
    return GrossVrabecQuatrupolar(components,pkg_param,model.references)
end


function a_polar(model,V,T,z,polar_model::GrossVrabecQuatrupolarModel,_data,_polar_data)
    μ̄² = model.params.dipole2.values
    Q̄² = model.params.quadrupole2.values
    has_dp = !all(iszero, μ̄²)
    has_qp = !all(iszero, Q̄²)
    _0 = zero(V+T+first(z))
    if !has_dp && !has_qp return _0 end
    a_mp_total = _0
    a_mp_total += has_dp && @f(a_dd,_data)
    a_mp_total += has_qp && @f(a_qq,_data)
    a_mp_total += has_dp && has_qp && @f(a_dq,_data)
    return a_mp_total
end

function a_dd(model ::GrossVrabecQuatrupolarModel, V, T, z, η=@f(data))
    a₂ = @f(a_2,η)
    iszero(a₂) && return zero(a₂)
    a₃ = @f(a_3,η)
    return a₂^2/(a₂-a₃)
end

function a_qq(model ::GrossVrabecQuatrupolarModel, V, T, z, η=@f(data))
    a₂ = @f(a_2_qq,η)
    iszero(a₂) && return zero(a₂)
    a₃ = @f(a_3_qq,η)
    return a₂^2/(a₂-a₃)
end

function a_dq(model ::GrossVrabecQuatrupolarModel, V, T, z, η=@f(data))
    a₂ = @f(a_2_dq,η)
    iszero(a₂) && return zero(a₂)
    a₃ = @f(a_3_dq,η)
    return a₂^2/(a₂-a₃)
end

# Deals only with DD & QQ
function a_2_qq(model ::QPCPSAFTModel, V, T, z, η=@f(data))
    Q̄² = model.params.quadrupole2.values
    ∑z = sum(z)
    ρ = N_A*∑z/V
    _a_2 = zero(T+V+first(z))
    m = model.params.segment.values
    ϵ = model.params.epsilon.values
    σ = model.params.sigma.values
    T⁻¹ = 1/T
    for i ∈ @comps
        ϵTii,zᵢ,Q̄²ᵢ = ϵ[i,i]*T⁻¹,z[i],Q̄²[i]
        iszero(Q̄²ᵢ) && continue
        _J2_ii = @f(J2,:QQ,i,i,η,m,ϵTii)
        _a_2 +=zᵢ^2*Q̄²ᵢ^2/σ[i,i]^7*_J2_ii
        for j ∈ 1:(i-1)
            Q̄²ⱼ = Q̄²[j]
            iszero(Q̄²ⱼ) && continue
            _J2_ij = @f(J2,:QQ,i,j,η,m,ϵ[i,j]*T⁻¹)
            _a_2 += 2*zᵢ*z[j]*Q̄²ᵢ*Q̄²[j]/σ[i,j]^7*_J2_ij
        end
    end
    _a_2 *= -π*9/16*ρ/(T*T)/(∑z*∑z)
    return _a_2
end

function a_2_dq(model ::QPCPSAFTModel, V, T, z, η=@f(data))
    
    Q̄² = model.params.quadrupole2.values
    μ̄² = model.params.dipole2.values
    ∑z = sum(z)
    ρ = N_A*∑z/V
    _a_2 = zero(T+V+first(z))
    m = model.params.segment.values
    ϵ = model.params.epsilon.values
    T⁻¹ = 1/T
    σ = model.params.sigma.values
    for i ∈ @comps
        ϵi,zi,μ̄²i,Q̄²i = ϵ[i,i],z[i],μ̄²[i],Q̄²[i]
        iszero(μ̄²i) && continue
        for j ∈ @comps
            ϵij,zj,μ̄²j,Q̄²j = ϵ[i,j],z[j],μ̄²[j],Q̄²[j]
            iszero(Q̄²j) && continue
            σij5 = σ[i,j]^5
            ϵ_TS = ϵij*T⁻¹
            _J2_ij = @f(J2,:DQ,i,j,η,m,ϵij*T⁻¹)
            _a_2 += zi*zj*μ̄²i*Q̄²j/σij5*_J2_ij
        end
    end
    _a_2 *= -π*9/4*ρ/(T*T)/(∑z*∑z)
    return _a_2
end

function a_3_qq(model ::QPCPSAFTModel, V, T, z, η=@f(data))
    _0 = zero(T+V+first(z))
    Q̄² = model.params.quadrupole2.values
    ∑z = sum(z)
    ρ = N_A*∑z/V
    _a_3 = _0
    m = model.params.segment.values
    σ = model.params.sigma.values
    nc = length(model)
    for i ∈ 1:nc
        zi,Q̄²i = z[i],Q̄²[i]
        iszero(Q̄²i) && continue
        _J3_iii = @f(J3,:QQ,i,i,i,η,m)
        a_3_i = zi*Q̄²i/σ[i,i]^3
        _a_3 += a_3_i^3*_J3_iii
        for j ∈ i+1:nc
            zj,Q̄²j = z[j],Q̄²[j]
            iszero(Q̄²j) && continue
            σij⁻³ = 1/σ[i,j]^3
            a_3_iij = zi*Q̄²i*σij⁻³
            a_3_ijj = zj*Q̄²j*σij⁻³
            a_3_j = zj*Q̄²j/σ[j,j]^3
            _J3_iij = @f(J3,:QQ,i,i,j,η,m)
            _J3_ijj = @f(J3,:QQ,i,j,j,η,m)
            _a_3 += 3*a_3_iij*a_3_ijj*(a_3_i*_J3_iij + a_3_j*_J3_ijj)
            for k ∈ j+1:nc
                zk,Q̄²k = z[k],Q̄²[k]
                iszero(Q̄²k) && continue
                _J3_ijk = @f(J3,:QQ,i,j,k,η,m)
                _a_3 += 6*zi*zj*zk*Q̄²i*Q̄²j*Q̄²k*σij⁻³/(σ[i,k]*σ[j,k])^3*_J3_ijk
            end
        end
    end
    _a_3 *= 9*π^2/16*ρ^2/(T*T*T)/(∑z*∑z*∑z)
    return _a_3
end

function a_3_dq(model ::QPCPSAFTModel, V, T, z, η=@f(data))
    μ̄² = model.params.dipole2.values
    Q̄² = model.params.quadrupole2.values
    ∑z = sum(z)
    ρ = N_A*∑z/V
    _a_3 = zero(T+V+first(z))
    m = model.params.segment.values
    σ = model.params.sigma.values
    nc = length(model)
    @inbounds for i ∈ 1:nc#dp_comps
        μ̄²i,zi = μ̄²[i],z[i]
        iszero(μ̄²i) && continue
        σi = σ[i,i]
        for j ∈ 1:nc #union(dq_comps,qp_comps)
            μ̄²j,zj,Q̄²j = μ̄²[j],z[j],Q̄²[j]
            iszero(μ̄²j) & iszero(Q̄²j) && continue
            σj = σ[j,j]
            σij = σ[i,j]
            for k ∈ 1:nc#qp_comps
                μ̄²k,zk,Q̄²k = μ̄²[k],z[k],Q̄²[k]
                iszero(Q̄²k) && continue
                _J3_ijk = @f(J3,:DQ,i,j,k,η,m)
                _a_3 += zi*zj*zk*σi/
                    (σ[k,k]*(σij*σ[i,k]*σ[j,k])^2)*
                    μ̄²i*Q̄²k*(σj*μ̄²j+1.19374/σj*Q̄²j)*_J3_ijk
            end
        end
    end
    _a_3 *= -ρ^2/(T*T*T)/(∑z*∑z*∑z)
    return _a_3
end

function J2(model::QPCPSAFTModel, V, T, z, type::Symbol, i, j, η = @f(data),m = model.params.segment.values,ϵT⁻¹ = model.params.epsilon.values[i,j]/T)
    m̄ = sqrt(m[i]*m[j])
    m̄1 = one(m̄)
    m̄2 = (m̄-1)/m̄
    m̄3 = m̄2*(m̄-2)/m̄
    m̄i = (m̄1,m̄2,m̄3)
    if type == :QQ
        aij_qq = ntuple(i -> dot(m̄i,QQ_consts.corr_a[i]),Val(5))
        bij_qq = ntuple(i -> dot(m̄i,QQ_consts.corr_b[i]),Val(5))
        #aij = NTuple{5,Float64}(m̄ for ai in QQ_consts.corr_a)
        #bij = NTuple{5,Float64}(dot(m̄i,bi) for bi in QQ_consts.corr_b)
        cij_qq = aij_qq .+ bij_qq .* ϵT⁻¹
        return evalpoly(η,cij_qq)
    else
        aij_dq = ntuple(i -> dot(m̄i,DQ_consts.corr_a[i]),Val(4))
        bij_dq = ntuple(i -> dot(m̄i,DQ_consts.corr_b[i]),Val(4))
        cij_dq = aij_dq .+ bij_dq .* ϵT⁻¹
        return evalpoly(η,cij_dq)
    end
end

function J3(model::QPCPSAFTModel, V, T, z, type::Symbol, i, j, k, η = @f(data),m = model.params.segment.values)
    m̄ = cbrt(m[i]*m[j]*m[k])
    m̄1 = one(m̄)
    m̄2 = (m̄-1)/m̄
    if type == :DQ
        m̄i_dq = (m̄1,m̄2)
        cijk_dq = ntuple(i -> dot(m̄i_dq,DQ_consts.corr_c[i]),Val(4))
        #cijk = NTuple{4}(dot(m̄i,ci) for ci in DQ_consts.corr_c)
        return evalpoly(η,cijk_dq)
    end
    m̄3 = m̄2*(m̄-2)/m̄
    m̄i = (m̄1,m̄2,m̄3)
    cijk_qq = ntuple(i -> dot(m̄i,QQ_consts.corr_c[i]),Val(4))
    return evalpoly(η,cijk_qq)
end

const QQ_consts = (
    corr_a =
    ((1.2378308, 1.2854109,	1.7942954),
    (2.4355031,	-11.465615,	0.7695103),
    (1.6330905,	22.086893,	7.2647923),
    (-1.6118152, 7.4691383,	94.486699),
    (6.9771185,	-17.197772,	-77.148458)),

    corr_b =
    ((0.4542718, -0.813734, 6.8682675),
    (-4.5016264, 10.06403, -5.1732238),
    (3.5858868,	-10.876631, -17.240207),
    (0., 0., 0.),
    (0., 0., 0.)),

    corr_c =
    ((-0.5000437, 2.0002094, 3.1358271),
    (6.5318692, -6.7838658, 7.2475888),
    (-16.01478, 20.383246, 3.0759478),
    (14.42597, -10.895984, 0.),
    (0., 0., 0.))
)

const DQ_consts = (
    corr_a =
    ((0.697095, -0.6734593, 0.6703408),
    (-0.6335541, -1.4258991, -4.3384718),
    (2.945509, 4.1944139, 7.2341684),
    (-1.4670273, 1.0266216, 0.)),

    corr_b =
    ((-0.4840383, 0.6765101, -1.1675601),
    (1.9704055, -3.0138675, 2.1348843),
    (-2.1185727, 0.4674266, 0.),
    (0., 0., 0.)),

    corr_c =
    ((0.7846431, -2.072202),
    (3.3427, -5.863904),
    (0.4689111, -0.1764887),
    (0., 0.))
)