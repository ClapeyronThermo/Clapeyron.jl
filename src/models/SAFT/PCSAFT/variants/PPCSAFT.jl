struct PPCSAFTParam <: EoSParam
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

abstract type PPCSAFTModel <: PCSAFTModel end
@newmodel PPCSAFT PPCSAFTModel PPCSAFTParam
default_references(::Type{PPCSAFT}) = ["10.1021/ie0003887", "10.1021/ie010954d"]
default_locations(::Type{PPCSAFT}) = ["SAFT/PCSAFT/PPCSAFT/","properties/molarmass.csv"]
function transform_params(::Type{PPCSAFT},params,components)
    sigma = params["sigma"]
    sigma.values .*= 1E-10
    params = saft_lorentz_berthelot(params)
    μ,m,Q = params["dipole"],params["segment"],params["quadrupole"]
    params["dipole2"] = SingleParam("Dipole squared",components, μ.^2 ./ m ./ k_B*1e-36*(1e-10*1e-3))
    params["quadrupole2"] = SingleParam("Quadrupole squared",components, Q.^2 ./ m ./ k_B*1e-56*(1e-10*1e-3))
    return params
end

"""
    PPCSAFTModel <: SAFTModel
    PPCSAFT(components;
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
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume
## Input models
- `idealmodel`: Ideal Model

## Description
Polar Perturbed-Chain SAFT (PPC-SAFT)

## References
1. Gross, J., & Vrabec, J. (2005). An equation-of-state contribution for polar components: Dipolar molecules. AIChE Journal, 52(3), 856-1282. [doi:10.1002/aic.10683](https://doi.org/10.1002/aic.10683)
2. Gross, J. (2005). An equation-of-state contribution for polar components: Quadrupolar molecules. AIChE Journal, 51(9), 2556-2568. [doi:10.1002/aic.10502](https://doi.org/10.1002/aic.10502)
3. Gross, J., & Vrabec, J. (2008). Vapor−Liquid Equilibria Simulation and an Equation of State Contribution for Dipole−Quadrupole Interactions. J. Phys. Chem. B, 112(1), 51-60. [doi:10.1021/jp072619u](https://doi.org/10.1021/jp072619u)
"""
PPCSAFT

export PPCSAFT

function recombine_impl!(model ::PPCSAFTModel)
    μ,Q,m = model.params.dipole,model.params.quadrupole,model.params.segment
    model.params.dipole2 .= μ.^2 ./ m ./ k_B * 1e-36*(1e-10*1e-3)
    model.params.quadrupole2 .= Q.^2 ./ m ./ k_B * 1e-56*(1e-10*1e-3)
    recombine_saft!(model)
end

function a_res(model ::PPCSAFTModel, V, T, z)
    _data = @f(data)
    return @f(a_hc,_data) + @f(a_disp,_data) + @f(a_assoc,_data) + @f(a_polar,_data)
end

function a_polar(model ::PPCSAFTModel, V, T, z, _data=@f(data))
    μ̄² = model.params.dipole2.values
    Q̄² = model.params.quadrupole2.values
    has_dp = !all(iszero, μ̄²)
    has_qp = !all(iszero, Q̄²)
    if !has_dp && !has_qp return 0. end
    a_p_total = 0.
    a_p_total += has_dp && @f(a_dd,_data)
    a_p_total += has_qp && @f(a_qq,_data)
    a_p_total += has_dp && has_qp && @f(a_dq,_data)
    return a_p_total
end

function a_mp(model ::PPCSAFTModel, V, T, z, type, _data=@f(data))
    a₂ = @f(a_2,type,_data)
    iszero(a₂) && return zero(a₂)
    a₃ = @f(a_3,type,_data)
    return a₂^2/(a₂-a₃)
end
function a_dd(model ::PPCSAFTModel, V, T, z, _data=@f(data))
    return @f(a_mp,"DD",_data)
end
function a_qq(model ::PPCSAFTModel, V, T, z, _data=@f(data))
    return @f(a_mp,"QQ",_data)
end
function a_dq(model ::PPCSAFTModel, V, T, z, _data=@f(data))
    a₂ = @f(a_2_dq,_data)
    iszero(a₂) && return zero(a₂)
    a₃ = @f(a_3_dq,_data)
    return a₂^2/(a₂-a₃)
end

# Deals only with DD & QQ
function a_2(model ::PPCSAFTModel, V, T, z, type, _data=@f(data))
    dp_comps, qp_comps = @f(polar_comps)
    P̄² = []
    p_comps = []
    p = 0
    coeff = 0.
    if type == "DD"
        if isempty(dp_comps) return 0. end
        P̄² = model.params.dipole2.values
        p_comps = dp_comps
        p = 3
        coeff = 1.
    end
    if type == "QQ"
        if isempty(qp_comps) return 0. end
        P̄² = model.params.quadrupole2.values
        p_comps = qp_comps
        p = 7
        coeff = 9/16
    end
    _,_,_,_,η,_ = _data
    ∑z = sum(z)
    ρ = N_A*∑z/V
    _a_2 = zero(T+V+first(z))
    # nc = length(model)
    m = model.params.segment.values
    ϵ = model.params.epsilon.values
    σ = model.params.sigma.values
    # μ̄² = model.params.dipole2.values
    # Q̄² = model.params.quadrupole2.values
    @inbounds for (idx, i) ∈ enumerate(p_comps)
        _J2_ii = @f(J2,type,i,i,η,m,ϵ)
        zᵢ = z[i]
        P̄²ᵢ = P̄²[i]
        _a_2 +=zᵢ^2*P̄²ᵢ^2/σ[i,i]^p*_J2_ii
        for j ∈ p_comps[idx+1:end]
            _J2_ij = @f(J2,type,i,j,η,m,ϵ)
            _a_2 += 2*zᵢ*z[j]*P̄²ᵢ*P̄²[j]/σ[i,j]^p*_J2_ij
        end
    end
    # @inbounds for i ∈ 1:nc
    #     _J2_ii = @f(J2,type,i,i,η,m,ϵ)
    #     zᵢ = z[i]
    #     μ̄²ᵢ = μ̄²[i]
    #     _a_2 +=zᵢ^2*μ̄²ᵢ^2/σ[i,i]^3*_J2_ii
    #     for j ∈ (i+1):nc
    #         _J2_ij = @f(J2,type,i,j,η,m,ϵ)
    #         _a_2 += 2*zᵢ*z[j]*μ̄²ᵢ*μ̄²[j]/σ[i,j]^3*_J2_ij
    #     end
    # end
    _a_2 *= -π*coeff*ρ/(T*T)/(∑z*∑z)
    return _a_2
end

function a_2_dq(model ::PPCSAFTModel, V, T, z, _data=@f(data))
    dp_comps, qp_comps = @f(polar_comps)
    Q̄² = model.params.quadrupole2.values
    μ̄² = model.params.dipole2.values
    _,_,_,_,η,_ = _data
    ∑z = sum(z)
    ρ = N_A*∑z/V
    _a_2 = zero(T+V+first(z))
    m = model.params.segment.values
    ϵ = model.params.epsilon.values
    σ = model.params.sigma.values
    @inbounds for i in dp_comps
        for j ∈ qp_comps
            _J2_ij = @f(J2,type,i,j,η,m,ϵ)
            _a_2 += z[i]*z[j]*μ̄²[i]*Q̄²[j]/σ[i,j]^5*_J2_ij
        end
    end
    _a_2 *= -π*9/4*ρ/(T*T)/(∑z*∑z)
    return _a_2
end

function a_3(model ::PPCSAFTModel, V, T, z, type, _data=@f(data))
    dp_comps, qp_comps = @f(polar_comps)
    P̄² = []
    p_comps = []
    p = 0
    coeff = 0.
    if type == "DD"
        if isempty(dp_comps) return 0. end
        P̄² = model.params.dipole2.values
        p_comps = dp_comps
        p = 1
        coeff = -4*π^2/3
    end
    if type == "QQ"
        if isempty(qp_comps) return 0. end
        P̄² = model.params.quadrupole2.values
        p_comps = qp_comps
        p = 3
        coeff = 9*π^2/16
    end

    ∑z = sum(z)
    ρ = N_A*∑z/V
    _,_,_,_,η,_ = _data
    _a_3 = zero(T+V+first(z))
    m = model.params.segment.values
    ϵ = model.params.epsilon.values
    σ = model.params.sigma.values
    # μ̄² = model.params.dipole2.values
    # nc = length(model)
    @inbounds for (idx_i,i) ∈ enumerate(p_comps)
        _J3_iii = @f(J3,type,i,i,i,η,m)
        zi,P̄²i = z[i],P̄²[i]
        a_3_i = zi*P̄²i/σ[i,i]^p
        _a_3 += a_3_i^3*_J3_iii
        for (idx_j,j) ∈ enumerate(p_comps[idx_i+1:end])
            zj,P̄²j = z[j],P̄²[j]
            σij⁻ᵖ = 1/σ[i,j]^p
            a_3_iij = zi*P̄²i*σij⁻ᵖ
            a_3_ijj = zj*P̄²j*σij⁻ᵖ
            a_3_j = zj*P̄²j/σ[j,j]^p
            _J3_iij = @f(J3,type,i,i,j,η,m)
            _J3_ijj = @f(J3,type,i,j,j,η,m)
            _a_3 += 3*a_3_iij*a_3_ijj*(a_3_i*_J3_iij + a_3_j*_J3_ijj)
            for k ∈ p_comps[idx_i+idx_j+1:end]
                zk,P̄²k = z[k],P̄²[k]
                _J3_ijk = @f(J3,type,i,j,k,η,m)
                _a_3 += 6*zi*zj*zk*P̄²i*P̄²j*P̄²k*σij⁻ᵖ/(σ[i,k]*σ[j,k])^p*_J3_ijk
            end
        end
    end
    # @inbounds for i ∈ 1:nc
    #     _J3_iii = @f(J3,type,i,i,i,η,m)
    #     zi,μ̄²i = z[i],μ̄²[i]
    #     a_3_i = zi*μ̄²i/σ[i,i]
    #     _a_3 += a_3_i^3*_J3_iii
    #     for j ∈ i+1:nc
    #         zj,μ̄²j = z[j],μ̄²[j]
    #         σij⁻¹ = 1/σ[i,j]
    #         a_3_iij = zi*μ̄²i*σij⁻¹
    #         a_3_ijj = zj*μ̄²j*σij⁻¹
    #         a_3_j = zj*μ̄²j/σ[j,j]
    #         _J3_iij = @f(J3,type,i,i,j,η,m)
    #         _J3_ijj = @f(J3,type,i,j,j,η,m)
    #         _a_3 += 3*a_3_iij*a_3_ijj*(a_3_i*_J3_iij + a_3_j*_J3_ijj)
    #         for k ∈ j+1:nc
    #             zk,μ̄²k = z[k],μ̄²[k]
    #             _J3_ijk = @f(J3,type,i,j,k,η,m)
    #             _a_3 += 6*zi*zj*zk*μ̄²i*μ̄²j*μ̄²k*σij⁻¹/σ[i,k]/σ[j,k]*_J3_ijk
    #         end
    #     end
    # end
    _a_3 *= coeff*ρ^2/(T*T*T)/(∑z*∑z*∑z)
    return _a_3
end


function a_3_original(model ::PPCSAFTModel, V, T, z, type, _data=@f(data))
    ∑z = sum(z)
    ρ = N_A*∑z/V
    _,_,_,_,η,_ = _data
    _a_3 = zero(T+V+first(z))
    m = model.params.segment.values
    ϵ = model.params.epsilon.values
    σ = model.params.sigma.values
    μ̄² = model.params.dipole2.values
    nc = length(model)

    @inbounds for i ∈ 1:nc
        _J3_iii = @f(J3,"DD",i,i,i,η,m)
        zi,μ̄²i = z[i],μ̄²[i]
        a_3_i = zi*μ̄²i/σ[i,i]
        _a_3 += a_3_i^3*_J3_iii
        for j ∈ i+1:nc
            zj,μ̄²j = z[j],μ̄²[j]
            σij⁻¹ = 1/σ[i,j]
            a_3_iij = zi*μ̄²i*σij⁻¹
            a_3_ijj = zj*μ̄²j*σij⁻¹
            a_3_j = zj*μ̄²j/σ[j,j]
            _J3_iij = @f(J3,"DD",i,i,j,η,m)
            _J3_ijj = @f(J3,"DD",i,j,j,η,m)
            _a_3 += 3*a_3_iij*a_3_ijj*(a_3_i*_J3_iij + a_3_j*_J3_ijj)
            for k ∈ j+1:nc
                zk,μ̄²k = z[k],μ̄²[k]
                _J3_ijk = @f(J3,"DD",i,j,k,η,m)
                _a_3 += 6*zi*zj*zk*μ̄²i*μ̄²j*μ̄²k*σij⁻¹/σ[i,k]/σ[j,k]*_J3_ijk
            end
        end
    end
    _a_3 *= -4π^2/3*ρ^2/(T*T*T)/(∑z*∑z*∑z)
    return _a_3
end

function a_3_dq(model ::PPCSAFTModel, V, T, z, _data=@f(data))
    dp_comps, qp_comps = @f(polar_comps)
    μ̄² = model.params.dipole2.values
    Q̄² = model.params.quadrupole2.values
    _,_,_,_,η,_ = _data
    ∑z = sum(z)
    ρ = N_A*∑z/V
    _a_3 = zero(T+V+first(z))
    m = model.params.segment.values
    ϵ = model.params.epsilon.values
    σ = model.params.sigma.values
    @inbounds for i ∈ dp_comps
        for j ∈ union(dp_comps, qp_comps)
            for k ∈ qp_comps
                _J3_ijk = @f(J3,type,i,j,k,η,m)
                _a_3 += z[i]*z[j]*z[k]*σ[i,i]/
                    (σ[k,k]*(σ[i,j]*σ[i,k]*σ[j,k])^2)*
                    μ̄²[i]*Q̄²[k]*(σ[j,j]*μ̄²[j]+1.19374/σ[j,j]*Q̄²[j])*_J3_ijk
            end
        end
    end
    _a_3 *= -ρ^2/(T*T*T)/(∑z*∑z*∑z)
    return _a_3
end

# Returns [Dipole Comp idxs], [Quadrupole Comp idxs]
function polar_comps(model, V, T, z)
    μ̄² = model.params.dipole2.values
    Q̄² = model.params.quadrupole2.values
    dipole_comps = []
    quadrupole_comps = []
    for i in @comps
        if !iszero(μ̄²[i]) push!(dipole_comps,i) end
        if !iszero(Q̄²[i]) push!(quadrupole_comps,i) end
    end
    return dipole_comps, quadrupole_comps
end


function J2(model, V, T, z, type, i, j, η = @f(ζ,3),m = model.params.segment.values,ϵ = model.params.epsilon.values)
    ϵT = ϵ[i,j]/T
    m̄ = sqrt(m[i]*m[j])
    if type == "DD" m̄ = min(m̄,2*one(η)) end
    consts = eval(Symbol(type*"_consts"))
    aij = _ppcsaft_corr_poly(consts.corr_a,m̄)
    bij = _ppcsaft_corr_poly(consts.corr_b,m̄)
    cij = aij .+ bij .* ϵT
    return evalpoly(η,cij)
end

function J3(model, V, T, z, type, i, j, k, η = @f(ζ,3),m = model.params.segment.values)
    corr_c = eval(Symbol(type*"_consts")).corr_c
    m̄ = cbrt(m[i]*m[j]*m[k])
    if type == "DD" m̄ = min(m̄,2*one(η)) end
    if type == "DQ"
        m̄1 = 1.
        m̄2 = (m̄-1)/m̄
        m̄i = (m̄1,m̄2)
        cijk = NTuple{5}(dot(m̄i,ci) for ci in corr)
        return evalpoly(η,cijk)
    end
    cijk = _ppcsaft_corr_poly(corr_c,m̄)
    return evalpoly(η,cijk)
end

function _ppcsaft_corr_poly(corr,m̄)
    m̄1 = 1.
    m̄2 = (m̄-1)/m̄
    m̄3 = m̄2*(m̄-2)/m̄
    m̄i = (m̄1,m̄2,m̄3)
    return NTuple{5}(dot(m̄i,ci) for ci in corr)
end

#Optimizations for Single Component PPCSAFT

function d(model::PPCSAFT, V, T, z::SingleComp)
    ϵ = only(model.params.epsilon.values)
    σ = only(model.params.sigma.values)
    return SA[σ*(1 - 0.12*exp(-3ϵ/T))]
end

const PPCSAFTconsts = (
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
    (0.,	0.,	0.)),

    corr_c =
    ((7.846431, -20.72202),
    (33.427, -58.63904),
    (4.689111, -1.764887),
    (0., 0.))
)