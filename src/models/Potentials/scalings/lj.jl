abstract type LJScalingModel <: EoSModel end

function σϵ_m_vdw1f(model::EoSModel,V,T,z)
    ϵ = model.params.epsilon.values
    σ = model.params.sigma.values
    m = model.params.segment.values
    return σϵ_m_vdw1f(ϵ,σ,m,V,T,z)
end

function σϵ_m_vdw1f(ϵ,σ,m,V,T,z)
    σϵ_mix = zero(eltype(z))
    σ_mix = zero(eltype(z))
    comps = length(z)
    ∑z = sum(z)
    m̄ = dot(m,z)/∑z
    m̄inv = 1/m̄
    for i ∈ 1:comps
        zi = z[i]*m[i]*m̄inv
        zii = zi*zi
        σ3 = σ[i,i]^3
        σϵ_mix += zii*σ3*ϵ[i,i]
        σ_mix += zii*σ3
        for j ∈ 1:(i-1)
            σ3ij = σ[i,j]^3
            zij = zi*z[j]*m[j]*m̄inv
            σϵ_mix += 2*zij*σ3ij*ϵ[i,j]
            σ_mix += 2*zij*σ3ij
        end
    end
    return σ_mix/(∑z*∑z), σϵ_mix/σ_mix, m̄
end

σϵ_m_vdw1f(ϵ,σ,V,T,z) = σϵ_m_vdw1f(ϵ,σ,FillArrays.Fill(one(eltype(σ)),length(z)),V,T,z)

struct LJScalingParam{T} <: ParametricEoSParam{T}
    Mw::SingleParam{T}
    sigma::PairParam{T}
    epsilon::PairParam{T}
end

function LJScalingParam(Mw,sigma,epsilon)
    return build_parametric_param(KolafaNezbedaParam,Mw,sigma,epsilon)
end

function transform_params_lj(params)
    sigma = params["sigma"]
    sigma.values .*= 1E-10
    return saft_lorentz_berthelot(params)
end

@newmodelsimple LJScaling LJScalingModel LJScalingParam{T}

function LJScaling(params::LJScalingParam{T}) where T
    comps = params.Mw.components
    return LJScaling(comps,params,comps) #use already allocated vector.
end

function shape_factors(scaling::LJScalingModel,ref,T,z)
    ϵ, σ = model.params.epsilon.values, model.params.sigma.values
    σ3, ϵ̄, _ = σϵ_m_vdw1f(ϵ,σ,1.0,T,z)
    f = ϵ̄
    h = N_A*σ3
    return f,h
end