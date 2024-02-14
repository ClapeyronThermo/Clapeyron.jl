struct PPCSAFTParam <: EoSParam
    Mw::SingleParam{Float64}
    segment::SingleParam{Float64}
    sigma::PairParam{Float64}
    epsilon::PairParam{Float64}
    dipole::SingleParam{Float64}
    dipole2::SingleParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
end

abstract type PPCSAFTModel <: PCSAFTModel end
@newmodel PPCSAFT PPCSAFTModel PPCSAFTParam

#=
at the moment, we use the Gross-Vrabec polar term
and the technically correct name is PCPSAFT
so, for the moment, we alias that.
=#
const PCPSAFTModel = PPCSAFTModel
const PCPSAFTParam = PPCSAFTParam
const PCPSAFT = PPCSAFT

default_references(::Type{PPCSAFT}) = ["10.1021/ie0003887", "10.1021/ie010954d"]
default_locations(::Type{PPCSAFT}) = ["SAFT/PCSAFT/PPCSAFT/","properties/molarmass.csv"]
function transform_params(::Type{PPCSAFT},params,components)
    sigma = params["sigma"]
    sigma.values .*= 1E-10
    params = saft_lorentz_berthelot(params)
    μ,m = params["dipole"],params["segment"]
    params["dipole2"] = SingleParam("Dipole squared",components, μ.^2 ./ m ./ k_B*1e-36*(1e-10*1e-3))
    return params
end

"""
    PPCSAFTModel <: PCSAFTModel

    const PCPSAFT = PPCSAFT
    PPCSAFT(components;
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
- `dipole`: Single Parameter (`Float64`) - Dipole moment `[D]`
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
"""
PPCSAFT
@doc (@doc PPCSAFT) PCPSAFT
export PPCSAFT,PCPSAFT

function recombine_impl!(model ::PPCSAFTModel)
    μ,m = model.params.dipole,model.params.segment
    model.params.dipole2 .= μ.^2 ./ m ./ k_B * 1e-36*(1e-10*1e-3)
    recombine_saft!(model)
end

function a_res(model ::PPCSAFTModel, V, T, z)
    _data = @f(data)
    return @f(a_hc,_data) + @f(a_disp,_data) + @f(a_assoc,_data) + @f(a_polar,_data)
end

function a_polar(model ::PPCSAFTModel, V, T, z, _data=@f(data))
    dipole = model.params.dipole.values
    if all(iszero,dipole)
        return zero(V+T+first(z))
    end
    a₂ = @f(a_2,_data)
    iszero(a₂) && return zero(a₂)
    a₃ = @f(a_3,_data)
    return a₂^2/(a₂-a₃)
end

function a_2(model ::PPCSAFTModel, V, T, z, _data=@f(data))
    _,_,_,_,η,_ = _data
    ∑z = sum(z)
    ρ = N_A*∑z/V
    _a_2 = zero(T+V+first(z))
    nc = length(model)
    m = model.params.segment.values
    ϵ = model.params.epsilon.values
    σ = model.params.sigma.values
    μ̄² = model.params.dipole2.values
    @inbounds for i ∈ 1:nc
        _J2_ii = @f(J2,i,i,η,m,ϵ)
        zᵢ = z[i]
        μ̄²ᵢ = μ̄²[i]
        iszero(μ̄²ᵢ) && continue
        _a_2 +=zᵢ^2*μ̄²ᵢ^2/σ[i,i]^3*_J2_ii
        for j ∈ (i+1):nc
            μ̄²ⱼ = μ̄²[j]
            iszero(μ̄²ⱼ) && continue
            _J2_ij = @f(J2,i,j,η,m,ϵ)
            _a_2 += 2*zᵢ*z[j]*μ̄²ᵢ*μ̄²ⱼ/σ[i,j]^3*_J2_ij
        end
    end
    _a_2 *= -π*ρ/(T*T)/(∑z*∑z)
    return _a_2
end

function a_3(model ::PPCSAFTModel, V, T, z, _data=@f(data))
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
        zi,μ̄²i = z[i],μ̄²[i]
        iszero(μ̄²i) && continue
        _J3_iii = @f(J3,i,i,i,η,m)
        a_3_i = zi*μ̄²i/σ[i,i]
        _a_3 += a_3_i^3*_J3_iii
        for j ∈ i+1:nc
            zj,μ̄²j = z[j],μ̄²[j]
            iszero(μ̄²j) && continue
            σij⁻¹ = 1/σ[i,j]
            a_3_iij = zi*μ̄²i*σij⁻¹
            a_3_ijj = zj*μ̄²j*σij⁻¹
            a_3_j = zj*μ̄²j/σ[j,j]
            _J3_iij = @f(J3,i,i,j,η,m)
            _J3_ijj = @f(J3,i,j,j,η,m)
            _a_3 += 3*a_3_iij*a_3_ijj*(a_3_i*_J3_iij + a_3_j*_J3_ijj)
            for k ∈ j+1:nc
                zk,μ̄²k = z[k],μ̄²[k]
                iszero(μ̄²k) && continue
                _J3_ijk = @f(J3,i,j,k,η,m)
                _a_3 += 6*zi*zj*zk*μ̄²i*μ̄²j*μ̄²k*σij⁻¹/σ[i,k]/σ[j,k]*_J3_ijk
            end
        end
    end
    _a_3 *= -4π^2/3*ρ^2/(T*T*T)/(∑z*∑z*∑z)
    return _a_3
end

function J2(model, V, T, z, i, j, η = @f(ζ,3),m = model.params.segment.values,ϵ = model.params.epsilon.values)
    corr_a = PPCSAFTconsts.corr_a
    corr_b = PPCSAFTconsts.corr_b
    m̄ = min(sqrt(m[i]*m[j]),2*one(m[i]))
    ϵT = ϵ[i,j]/T
    aij = _ppcsaft_corr_poly(corr_a,m̄)
    bij = _ppcsaft_corr_poly(corr_b,m̄)
    cij = aij .+ bij .* ϵT
    return evalpoly(η,cij)
end

function J3(model, V, T, z, i, j, k, η = @f(ζ,3),m = model.params.segment.values)
    corr_c = PPCSAFTconsts.corr_c
    m̄ = min(cbrt(m[i]*m[j]*m[k]),2*one(m[i]))
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