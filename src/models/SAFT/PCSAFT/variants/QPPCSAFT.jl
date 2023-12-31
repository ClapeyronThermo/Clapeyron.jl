struct QPPCSAFTParam <: EoSParam
    Mw::SingleParam{Float64}
    segment::SingleParam{Float64}
    sigma::PairParam{Float64}
    epsilon::PairParam{Float64}
    dipole::SingleParam{Float64}
    dipole2::SingleParam{Float64}
    quadrupole::SingleParam{Float64}
    quadrupole2::SingleParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
end

abstract type QPPCSAFTModel <: PPCSAFTModel end
@newmodel QPPCSAFT QPPCSAFTModel QPPCSAFTParam
default_references(::Type{QPPCSAFT}) = ["10.1002/aic.10502","10.1021/jp072619u"]
default_locations(::Type{QPPCSAFT}) = ["SAFT/PCSAFT/QPPCSAFT/","properties/molarmass.csv"] # Needs to add data for QPPCSAFT
function transform_params(::Type{QPPCSAFT},params,components)
    sigma = params["sigma"]
    sigma.values .*= 1E-10
    params = saft_lorentz_berthelot(params)
    μ,m,Q = params["dipole"],params["segment"],params["quadrupole"]
    params["dipole2"] = SingleParam("Dipole squared",components, μ.^2 ./ m ./ k_B*1e-36*(1e-10*1e-3))
    params["quadrupole2"] = SingleParam("Quadrupole squared",components, Q.^2 ./ m ./ k_B*1e-56*(1e-10*1e-3))
    return params
end

"""
    QPPCSAFTModel <: SAFTModel
    QPPCSAFT(components;
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
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Paramater (no units)
- `dipole`: Single Parameter (`Float64`) - Dipole moment `[D]`
- `quadrupole`: Single Parameter (`Float64`) - Quadrupole moment `[DA°]`
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m^3]`

## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Pair Parameter (`Float64`) - Mixed segment Diameter `[m]`
- `epsilon`: Pair Parameter (`Float64`) - Mixed reduced dispersion energy`[K]`
- `dipole`: Single Parameter (`Float64`) - Dipole moment `[D]`
- `quadrupole`: Single Parameter (`Float64`) - Quadrupole moment `[DA°]`
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume
## Input models
- `idealmodel`: Ideal Model

## Description
Polar Perturbed-Chain SAFT, including Quadrupolar interaction contributions (QPPC-SAFT)

## References
1. Gross, J. (2005). An equation-of-state contribution for polar components: Quadrupolar molecules. AIChE Journal, 51(9), 2556-2568. [doi:10.1002/aic.10502](https://doi.org/10.1002/aic.10502)
2. Gross, J., & Vrabec, J. (2008). Vapor−Liquid Equilibria Simulation and an Equation of State Contribution for Dipole−Quadrupole Interactions. J. Phys. Chem. B, 112(1), 51-60. [doi:10.1021/jp072619u](https://doi.org/10.1021/jp072619u)
"""
QPPCSAFT

export QPPCSAFT

function recombine_impl!(model ::QPPCSAFTModel)
    μ,Q,m = model.params.dipole,model.params.quadrupole,model.params.segment
    model.params.dipole2 .= μ.^2 ./ m ./ k_B * 1e-36*(1e-10*1e-3)
    model.params.quadrupole2 .= Q.^2 ./ m ./ k_B * 1e-56*(1e-10*1e-3)
    recombine_saft!(model)
end

function a_res(model ::QPPCSAFTModel, V, T, z)
    _data = @f(data)
    return @f(a_hc,_data) + @f(a_disp,_data) + @f(a_assoc,_data) + @f(a_mp,_data)
end

function a_mp(model ::QPPCSAFTModel, V, T, z, _data=@f(data))
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

function a_dd(model ::QPPCSAFTModel, V, T, z, _data=@f(data))
    return @f(a_polar,_data)
end

function a_qq(model ::QPPCSAFTModel, V, T, z, _data=@f(data))
    a₂ = @f(a_2_qq,_data)
    iszero(a₂) && return zero(a₂)
    a₃ = @f(a_3_qq,_data)
    return a₂^2/(a₂-a₃)
end
function a_dq(model ::QPPCSAFTModel, V, T, z, _data=@f(data))
    a₂ = @f(a_2_dq,_data)
    iszero(a₂) && return zero(a₂)
    a₃ = @f(a_3_dq,_data)
    return a₂^2/(a₂-a₃)
end

# Deals only with DD & QQ
function a_2_qq(model ::QPPCSAFTModel, V, T, z, _data=@f(data))
    Q̄² = model.params.quadrupole2.values
    _,_,_,_,η,_ = _data
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

function a_2_dq(model ::QPPCSAFTModel, V, T, z, _data=@f(data))
    
    Q̄² = model.params.quadrupole2.values
    μ̄² = model.params.dipole2.values
    _,_,_,_,η,_ = _data
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
            ϵ_TS =  ϵij*T⁻¹
            _J2_ij = @f(J2,:DQ,i,j,η,m,ϵij*T⁻¹)
            _a_2 += zi*zj*μ̄²i*Q̄²j/σij5*_J2_ij
        end
    end
    _a_2 *= -π*9/4*ρ/(T*T)/(∑z*∑z)
    return _a_2
    dp_comps, qp_comps = @f(polar_comps)
    Q̄² = model.params.quadrupole2.values
    μ̄² = model.params.dipole2.values
    _,_,_,_,η,_ = _data
    ∑z = sum(z)
    ρ = N_A*∑z/V
    _a_2 = zero(T+V+first(z))
    m = model.params.segment.values
    ϵ = model.params.epsilon.values
    #ϵ_TS = [sqrt(ϵ[i,i]*ϵ[j,j]) for i ∈ @comps, j ∈ @comps]
    σ = model.params.sigma.values
    @inbounds for i ∈ dp_comps
        for j ∈ qp_comps
            _J2_ij = @f(J2,:DQ,i,j,η,m)
            _a_2 += z[i]*z[j]*μ̄²[i]*Q̄²[j]/σ[i,j]^5*_J2_ij
        end
    end
    _a_2 *= -π*9/4*ρ/(T*T)/(∑z*∑z)
    return _a_2
end

function a_3_qq(model ::QPPCSAFTModel, V, T, z, _data=@f(data))
    _0 = zero(T+V+first(z))
    Q̄² = model.params.quadrupole2.values
    ∑z = sum(z)
    ρ = N_A*∑z/V
    _,_,_,_,η,_ = _data
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

function a_3_dq(model ::QPPCSAFTModel, V, T, z, _data=@f(data))
    μ̄² = model.params.dipole2.values
    Q̄² = model.params.quadrupole2.values
    _,_,_,_,η,_ = _data
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

function polar_comps(model, V, T, z)
    μ̄² = model.params.dipole2.values
    Q̄² = model.params.quadrupole2.values
    dipole_comps = Int[]
    quadrupole_comps = Int[]
    for i ∈ @comps
        if !iszero(μ̄²[i]) push!(dipole_comps,i) end
        if !iszero(Q̄²[i]) push!(quadrupole_comps,i) end
    end
    return dipole_comps, quadrupole_comps
end


function J2(model::QPPCSAFTModel, V, T, z, type::Symbol, i, j, η = @f(ζ,3),m = model.params.segment.values,ϵT⁻¹ = model.params.epsilon.values[i,j]/T)
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

function J3(model::QPPCSAFTModel, V, T, z, type::Symbol, i, j, k, η = @f(ζ,3),m = model.params.segment.values)
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