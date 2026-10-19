struct ChainLJScalingParam{T} <: ParametricEoSParam{T}
    Mw::SingleParam{T}
    segment::SingleParam{T}
    sigma::PairParam{T}
    epsilon::PairParam{T}
end

function ChainLJScalingParam(Mw,segment,sigma,epsilon)
    return build_parametric_param(ChainLJScalingParam,Mw,segment,sigma,epsilon)
end

@newmodelsimple ChainLJScaling LJScalingModel LJScalingParam{T}

function ChainLJScaling(params::ChainLJScalingParam{T}) where T
    comps = params.Mw.components
    return ChainLJScaling(comps,params,comps) #use already allocated vector.
end

function shape_factors(scaling::ChainLJScaling,ref,T,z)
    ϵ, σ, m = model.params.epsilon.values, model.params.sigma.values, model.params.segment.values
    σ3, ϵ̄, m̄ = σϵ_m_vdw1f(ϵ,σ,m,1.0,T,z)
    f = ϵ̄
    h = m̄*N_A*σ3
    return f,h
end