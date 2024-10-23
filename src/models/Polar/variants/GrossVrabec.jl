#=
abstract type GrossVrabecModel <: PolarModel end

struct GrossVrabecParam{T} <: ParametricEoSParam{T}
    sigma::PairParam{T}
    epsilon::PairParam{T}
    segment::SingleParam{T}
    dipole::SingleParam{T}
    dipole2::SingleParam{T}
end

@newmodelsimple GrossVrabec GrossVrabecModel GrossVrabecParam

function GrossVrabec(model::EoSModel)
    components = model.components
    sigma = model.params.sigma
    epsilon = model.params.epsilon
    segment = model.params.segment
    dipole = model.params.dipole
    paramtype = typeof(model.params)
    if hasfield(paramtype,:dipole2)
        dipole2 = model.params.dipole
    else
        dipole2 = μ.^2 ./ m ./ k_B*1e-36*(1e-10*1e-3)
    end
    pkg_param = GrossVrabecParam(sigma,epsilon,segment,dipole,dipole2)
    return GrossVrabec(components,pkg_param,model.references)
end

function d(model::GrossVrabec,V,T,z)
    return diagvalues(model.sigma.values)
end

function polar_data(model,V,T,z,polar_model::GrossVrabecModel,_data)
    eta0 = packing_fraction(model,_data) #try getter
    if isnothing(eta0)
        return @f(packing_fraction) #calculate packing fraction if not available.
    end
end

function a_polar(model,V,T,z,polar_model::GrossVrabecModel,_data,_polar_data)
    dipole = model.params.dipole.values
    if all(iszero,dipole)
        return zero(V+T+first(z))
    end
    a₂ = @f(a_2_dd,_polar_data)
    iszero(a₂) && return zero(a₂)
    a₃ = @f(a_2_dd,_polar_data)
    return a₂^2/(a₂-a₃)
end

function a_2_dd(model::GrossVrabecModel, V, T, z, η = @f(packing_fraction))
    ∑z = sum(z)
    ρ = N_A*∑z/V
    _a_2 = zero(T+V+first(z))
    m = model.params.segment.values
    ϵ = model.params.sigma.values
    σ = model.params.epsilon.values
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

function a_3_dd(model::GrossVrabecModel, V, T, z, η = @f(packing_fraction))
    ∑z = sum(z)
    ρ = N_A*∑z/V
    _a_3 = zero(T+V+first(z))
    m = model.params.segment.values
    ϵ = model.params.sigma.values
    σ = model.params.epsilon.values
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
            _J3_iij = @f(J3,i,i,j,η,m,ϵ)
            _J3_ijj = @f(J3,i,j,j,η,m,ϵ)
            _a_3 += 3*a_3_iij*a_3_ijj*(a_3_i*_J3_iij + a_3_j*_J3_ijj)
            for k ∈ j+1:nc
                zk,μ̄²k = z[k],μ̄²[k]
                iszero(μ̄²k) && continue
                _J3_ijk = @f(J3,i,j,k,η,m,ϵ)
                _a_3 += 6*zi*zj*zk*μ̄²i*μ̄²j*μ̄²k*σij⁻¹/σ[i,k]/σ[j,k]*_J3_ijk
            end
        end
    end
    _a_3 *= -4π^2/3*ρ^2/(T*T*T)/(∑z*∑z*∑z)
    return _a_3
end

function J2(model::GrossVrabecModel, V, T, z, i, j, η = @f(packing_fraction), m = model.params.segment.values,ϵ = model.params.epsilon.values)
    corr_a = GrossVrabecConsts.corr_a
    corr_b = GrossVrabecConsts.corr_b
    m̄ = min(sqrt(m[i]*m[j]),2*one(m[i]))
    ϵT = ϵ[i,j]/T
    aij = gv_corr_poly(corr_a,m̄)
    bij = gv_corr_poly(corr_b,m̄)
    cij = aij .+ bij .* ϵT
    return evalpoly(η,cij)
end

function J3(model::GrossVrabecModel, V, T, z, i, j, k, η = @f(packing_fraction), m = model.segment.values,ϵ = model.params.epsilon.values)
    corr_c = GrossVrabecConsts.corr_c
    m̄ = min(cbrt(m[i]*m[j]*m[k]),2*one(m[i]))
    cijk = gv_corr_poly(corr_c,m̄)
    return evalpoly(η,cijk)
end

function gv_corr_poly(corr,m̄)
    m̄1 = 1.
    m̄2 = (m̄-1)/m̄
    m̄3 = m̄2*(m̄-2)/m̄
    m̄i = (m̄1,m̄2,m̄3)
    return NTuple{5}(dot(m̄i,ci) for ci in corr)
end

const GrossVrabecConsts = (
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
=#
