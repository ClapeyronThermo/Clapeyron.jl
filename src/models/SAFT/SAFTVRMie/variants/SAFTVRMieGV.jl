"""
SAFT-VR Mie multipolar approach using GV polar terms.

Polar terms copied verbatim from PCPSAFT and QPCPSAFT since they are the same ones used here.
"""
struct SAFTVRMieGVParam{T} <: ParametricEoSParam{T}
    Mw::SingleParam{Float64}
    segment::SingleParam{Float64}
    sigma::PairParam{Float64}
    lambda_a::PairParam{Float64}
    lambda_r::PairParam{Float64}
    epsilon::PairParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
    dipole::SingleParam{Float64}
    dipole2::SingleParam{Float64}
    np::SingleParam{Float64}
    quadrupole::SingleParam{Float64}
    quadrupole2::SingleParam{Float64}
	nQ::SingleParam{Float64}   
end


function SAFTVRMieGVParam(Mw,segment,sigma,lambda_a,lambda_r,epsilon,epsilon_assoc,bondvol,dipole,dipole2,np,quadrupole,quadrupole2,nQ)
    return build_parametric_param(SAFTVRMieGVParam,Mw,segment,sigma,lambda_a,lambda_r,epsilon,epsilon_assoc,bondvol,dipole,dipole2,np,quadrupole,quadrupole2,nQ) 
end

abstract type SAFTVRMieGVModel <: SAFTVRMieModel end
@newmodel SAFTVRMieGV SAFTVRMieGVModel SAFTVRMieGVParam{T}
default_references(::Type{SAFTVRMieGV}) = ["10.1016/j.fluid.2017.09.027","10.1021/acs.jced.0c00705"]
default_locations(::Type{SAFTVRMieGV}) = ["SAFT/SAFTVRMie/SAFTVRMieGV/","properties/molarmass.csv"] 

function transform_params(::Type{SAFTVRMieGV},params,components)
    sigma = params["sigma"]
    sigma.values .*= 1E-10
    sigma = sigma_LorentzBerthelot(sigma)
    epsilon = epsilon_HudsenMcCoubreysqrt(params["epsilon"], sigma)
    lambda_a = lambda_LorentzBerthelot(params["lambda_a"])
    lambda_r = lambda_LorentzBerthelot(params["lambda_r"])
    m = params["segment"]

    μ = get!(params,"dipole") do
        SingleParam("dipole",components)
    end
    np = get!(params,"np") do
        SingleParam("np",components)
    end
    Q = get!(params,"quadrupole") do
        SingleParam("quadrupole",components)
    end
    nQ = get!(params,"nQ") do
        SingleParam("nQ",components)
    end
    
    params["sigma"] = sigma
    params["epsilon"] = epsilon
    params["lambda_a"] = lambda_a
    params["lambda_r"] = lambda_r
    params["dipole2"] = SingleParam("Dipole squared",components, np .* μ.^2 ./ m ./ k_B*1e-36*(1e-10*1e-3))
    params["quadrupole2"] = SingleParam("Quadrupole squared",components, nQ .* Q.^2 ./ m ./ k_B*1e-56*(1e-10*1e-3))
    return params
end

"""
    SAFTVRMieGVModel <: SAFTModel
    SAFTVRMieGV(components;
    idealmodel=BasicIdeal,
    userlocations=String[],
    ideal_userlocations=String[],
    verbose=false,
    assoc_options = AssocOptions()) # SS: is it possible to add another parameter here to select/deselect DQ term?

## Input parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Single Parameter (`Float64`) - Segment Diameter `[Å]`
- `epsilon`: Single Parameter (`Float64`) - Reduced dispersion energy `[K]`
- `lambda_a`: Pair Parameter (`Float64`) - Attractive range parameter (no units)
- `lambda_r`: Pair Parameter (`Float64`) - Repulsive range parameter (no units)
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Parameter (no units)
- `dipole`: Single Parameter (`Float64`) - Dipole moment `[D]`
- `np` : Single Parameter (`Float64`) - number of dipolar segments (no units)
- `quadrupole`: Single Parameter (`Float64`) - Quadrupole moment `[D·Å]`
- `nQ` : Single Parameter (`Float64`) - number of quadrupolar segments (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m³]`

## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Pair Parameter (`Float64`) - Mixed segment Diameter `[m]`
- `epsilon`: Pair Parameter (`Float64`) - Mixed reduced dispersion energy`[K]`
- `lambda_a`: Pair Parameter (`Float64`) - Attractive range parameter (no units)
- `lambda_r`: Pair Parameter (`Float64`) - Repulsive range parameter (no units)
- `dipole`: Single Parameter (`Float64`) - Dipole moment `[D]`
- `np` : Single Parameter (`Float64`) - number of dipolar segments (no units)
- `quadrupole`: Single Parameter (`Float64`) - Quadrupole moment `[D·Å]`
- `nQ` : Single Parameter (`Float64`) - number of quadrupolar segments (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m³]`

## Input models
- `idealmodel`: Ideal Model

## Description
Polar SAFT-VR with Mie potential, including Dipolar and Quadrupolar interaction contributions (SAFT-VR Mie GV)

## References
Polar terms:
1. Gross, J. (2005). An equation-of-state contribution for polar components: Quadrupolar molecules. AIChE Journal, 51(9), 2556-2568. [doi:10.1002/aic.10502](https://doi.org/10.1002/aic.10502)
2. Gross, J., & Vrabec, J. (2006). An equation-of-state contribution for polar components: Dipolar molecules. AIChE Journal, 52(3), 856-1282. [doi:10.1002/aic.10683](https://doi.org/10.1002/aic.10683)
3. Gross, J., & Vrabec, J. (2008). Vapor−Liquid Equilibria Simulation and an Equation of State Contribution for Dipole−Quadrupole Interactions. J. Phys. Chem. B, 112(1), 51-60. [doi:10.1021/jp072619u](https://doi.org/10.1021/jp072619u)

Implemented in SAFT-VR Mie:
1. Cripwell et al. (2017). SAFT-VR-Mie with an incorporated polar term for accurate holistic prediction of the thermodynamic properties of polar components. Fluid Phase Equilibria, 455, 24-42. [doi:10.1016/j.fluid.2017.09.027](https://doi.org/10.1016/j.fluid.2017.09.027)
2. Smith et al. (2020). A Quadrupolar SAFT-VR Mie Approach to Modeling Binary Mixtures of CO2 or Benzene with n-Alkanes or 1-Alkanols. J. Chem. Eng. Data, 65(12), 5778-5800. [doi:10.1021/acs.jced.0c00705](https://doi.org/10.1021/acs.jced.0c00705)
"""
SAFTVRMieGV

export SAFTVRMieGV

# SS: used to recalculate mix params if one param has changed. Copied verbatim from SAFTVRMie.jl. 
function recombine_impl!(model ::SAFTVRMieGVModel)
    μ,Q,m,np,nQ = model.params.dipole,model.params.quadrupole,model.params.segment,model.params.np,model.params.nQ
    model.params.dipole2 .= np .* μ.^2 ./ m ./ k_B * 1e-36*(1e-10*1e-3)  # 1e-49
    model.params.quadrupole2 .= nQ .* Q.^2 ./ m ./ k_B * 1e-56*(1e-10*1e-3)  # 1e-69
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

function a_res(model ::SAFTVRMieGVModel, V, T, z) 
    _data = @f(data) 
    return @f(a_hs,_data) + @f(a_dispchain,_data) + @f(a_assoc,_data) + @f(a_mp,_data)
end

function a_mp(model ::SAFTVRMieGVModel, V, T, z, _data=@f(data))
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

function a_dd(model ::SAFTVRMieGVModel, V, T, z, _data=@f(data))
    a₂ = @f(a_2_dd,_data)
    iszero(a₂) && return zero(a₂)
    a₃ = @f(a_3_dd,_data)
    return a₂^2/(a₂-a₃)
end
function a_qq(model ::SAFTVRMieGVModel, V, T, z, _data=@f(data))
    a₂ = @f(a_2_qq,_data)
    iszero(a₂) && return zero(a₂)
    a₃ = @f(a_3_qq,_data)
    return a₂^2/(a₂-a₃)
end
function a_dq(model ::SAFTVRMieGVModel, V, T, z, _data=@f(data))
    a₂ = @f(a_2_dq,_data)
    iszero(a₂) && return zero(a₂)
    a₃ = @f(a_3_dq,_data)
    return a₂^2/(a₂-a₃)
end

# Dipolar term
function a_2_dd(model ::SAFTVRMieGVModel, V, T, z, _data=@f(data))
    _,_,ζi,_,_,_,_ = _data
    _,_,_,η = ζi

    μ̄² = model.params.dipole2.values
    ∑z = sum(z)
    ρ = N_A*∑z/V
    _a_2 = zero(T+V+first(z))
    m = model.params.segment.values
    ϵ = model.params.epsilon.values
    σ = model.params.sigma.values
    T⁻¹ = 1/T
    for i ∈ @comps
        ϵTii,zᵢ,μ̄²ᵢ = ϵ[i,i]*T⁻¹,z[i],μ̄²[i]
        iszero(μ̄²ᵢ) && continue
        _J2_ii = @f(J2,:DD,i,i,η,m,ϵTii)
        _a_2 +=zᵢ^2*μ̄²ᵢ^2/σ[i,i]^3*_J2_ii
        for j ∈ 1:(i-1)
            μ̄²ⱼ = μ̄²[j]
            iszero(μ̄²ⱼ) && continue
            _J2_ij = @f(J2,:DD,i,j,η,m,ϵ[i,j]*T⁻¹)
            _a_2 += 2*zᵢ*z[j]*μ̄²ᵢ*μ̄²ⱼ/σ[i,j]^3*_J2_ij
        end
    end
    _a_2 *= -π*ρ/(T*T)/(∑z*∑z)
    return _a_2
end

function a_3_dd(model ::SAFTVRMieGVModel, V, T, z, _data=@f(data))
    _,_,ζi,_,_,_,_ = _data
    _,_,_,η = ζi

    _0 = zero(T+V+first(z))
    μ̄² = model.params.dipole2.values
    ∑z = sum(z)
    ρ = N_A*∑z/V
    _a_3 = _0
    m = model.params.segment.values
    σ = model.params.sigma.values
    nc = length(model)

    for i ∈ 1:nc
        zi,μ̄²i = z[i],μ̄²[i]
        iszero(μ̄²i) && continue
        _J3_iii = @f(J3,:DD,i,i,i,η,m)
        a_3_i = zi*μ̄²i/σ[i,i]
        _a_3 += a_3_i^3*_J3_iii
        for j ∈ i+1:nc
            zj,μ̄²j = z[j],μ̄²[j]
            iszero(μ̄²j) && continue
            σij⁻¹ = 1/σ[i,j]
            a_3_iij = zi*μ̄²i*σij⁻¹
            a_3_ijj = zj*μ̄²j*σij⁻¹
            a_3_j = zj*μ̄²j/σ[j,j]
            _J3_iij = @f(J3,:DD,i,i,j,η,m)
            _J3_ijj = @f(J3,:DD,i,j,j,η,m)
            _a_3 += 3*a_3_iij*a_3_ijj*(a_3_i*_J3_iij + a_3_j*_J3_ijj)
            for k ∈ j+1:nc
                zk,μ̄²k = z[k],μ̄²[k]
                iszero(μ̄²k) && continue
                _J3_ijk = @f(J3,:DD,i,j,k,η,m)
                _a_3 += 6*zi*zj*zk*μ̄²i*μ̄²j*μ̄²k*σij⁻¹/σ[i,k]/σ[j,k]*_J3_ijk
            end
        end
    end
    _a_3 *= -4π^2/3*ρ^2/(T*T*T)/(∑z*∑z*∑z)
    return _a_3
end


# Quadrupolar term
function a_2_qq(model ::SAFTVRMieGVModel, V, T, z, _data=@f(data))
    _,_,ζi,_,_,_,_ = _data
    _,_,_,η = ζi

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

function a_3_qq(model ::SAFTVRMieGVModel, V, T, z, _data=@f(data))
    _,_,ζi,_,_,_,_ = _data
    _,_,_,η = ζi

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

# Cross dipole-quadrupole term
function a_2_dq(model ::SAFTVRMieGVModel, V, T, z, _data=@f(data))
    _,_,ζi,_,_,_,_ = _data
    _,_,_,η = ζi
    
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
            ϵ_TS =  ϵij*T⁻¹
            _J2_ij = @f(J2,:DQ,i,j,η,m,ϵij*T⁻¹)
            _a_2 += zi*zj*μ̄²i*Q̄²j/σij5*_J2_ij
        end
    end
    _a_2 *= -π*9/4*ρ/(T*T)/(∑z*∑z)
    return _a_2
end

function a_3_dq(model ::SAFTVRMieGVModel, V, T, z, _data=@f(data))
    _,_,ζi,_,_,_,_ = _data
    _,_,_,η = ζi

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
            for k ∈ 1:nc #qp_comps
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
#=
function polar_comp(model, V, T, z)
    μ̄² = model.params.dipole2.values
    Q̄² = model.params.quadrupole2.values
    dipole_comps = Int[]
    quadrupole_comps = Int[]
    for i ∈ @comps
        if !iszero(μ̄²[i]) push!(dipole_comps,i) end
        if !iszero(Q̄²[i]) push!(quadrupole_comps,i) end
    end
    return dipole_comps, quadrupole_comps
end=#

function J2(model::SAFTVRMieGVModel, V, T, z, type::Symbol, i, j, η = @f(ζ0123,4), m = model.params.segment.values,ϵT⁻¹ = model.params.epsilon.values[i,j]/T)
    m̄ = sqrt(m[i]*m[j])
    m̄1 = one(m̄)
    m̄2 = (m̄-1)/m̄
    m̄3 = m̄2*(m̄-2)/m̄
    m̄i = (m̄1,m̄2,m̄3)
    if type == :QQ
        aij_qq = ntuple(i -> dot(m̄i,QQconsts.corr_a[i]),Val(5))
        bij_qq = ntuple(i -> dot(m̄i,QQconsts.corr_b[i]),Val(5))
        cij_qq = aij_qq .+ bij_qq .* ϵT⁻¹
        return evalpoly(η,cij_qq)
    elseif type == :DD
        # For DD, there is a restriction on the value of m̄
        m̄ = min(sqrt(m[i]*m[j]),2*one(m[i]))
        m̄1 = one(m̄)
        m̄2 = (m̄-1)/m̄
        m̄3 = m̄2*(m̄-2)/m̄
        m̄i = (m̄1,m̄2,m̄3)

        aij_dd = ntuple(i -> dot(m̄i,DD_consts.corr_a[i]),Val(5))
        bij_dd = ntuple(i -> dot(m̄i,DD_consts.corr_b[i]),Val(5))
        cij_dd = aij_dd .+ bij_dd .* ϵT⁻¹
        return evalpoly(η,cij_dd)
    else
        aij_dq = ntuple(i -> dot(m̄i,DQ_consts.corr_a[i]),Val(4))
        bij_dq = ntuple(i -> dot(m̄i,DQ_consts.corr_b[i]),Val(4))
        cij_dq = aij_dq .+ bij_dq .* ϵT⁻¹
        return evalpoly(η,cij_dq)
    end
end

function J3(model::SAFTVRMieGVModel, V, T, z, type::Symbol, i, j, k, η = @f(ζ0123,4), m = model.params.segment.values)
    m̄ = cbrt(m[i]*m[j]*m[k])
    m̄1 = one(m̄)
    m̄2 = (m̄-1)/m̄
    if type == :DQ
        m̄i_dq = (m̄1,m̄2)
        cijk_dq = ntuple(i -> dot(m̄i_dq,DQ_consts.corr_c[i]),Val(4))
        return evalpoly(η,cijk_dq)
    elseif type == :QQ
        m̄3 = m̄2*(m̄-2)/m̄
        m̄i = (m̄1,m̄2,m̄3)
        cijk_qq = ntuple(i -> dot(m̄i,QQ_consts.corr_c[i]),Val(4))
        return evalpoly(η,cijk_qq)
    else
        # For DD, there is a restriction on the value of m̄
        m̄ = min(cbrt(m[i]*m[j]*m[k]),2*one(m[i]))
        m̄1 = one(m̄)
        m̄2 = (m̄-1)/m̄
        m̄3 = m̄2*(m̄-2)/m̄
        m̄i = (m̄1,m̄2,m̄3)
        
        cijk_dd = ntuple(i -> dot(m̄i,DD_consts.corr_c[i]),Val(4))
        return evalpoly(η,cijk_dd)
    end

end


const DD_consts = (
    corr_a =
    ((0.3043504,0.9534641,-1.161008),
    (-0.1358588,-1.8396383,4.5258607),
    (1.4493329,2.013118,0.9751222),
    (0.3556977,-7.3724958,-12.281038),
    (-2.0653308,8.2374135,5.9397575)),

    corr_b =
    ((0.2187939,-0.5873164,3.4869576),
    (-1.1896431,1.2489132,-14.915974),
    (1.1626889,-0.508528,15.372022),
    (0.,0.,0.),
    (0.,0.,0.)),

    corr_c =
    ((-0.0646774,-0.9520876,-0.6260979),
    (0.1975882,2.9924258,1.2924686),
    (-0.8087562,-2.3802636,1.6542783),
    (0.6902849,-0.2701261,-3.4396744),
    (0.,0.,0.))
)

const QQconsts = (
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