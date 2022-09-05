abstract type CubicPlusChainModel <: EoSModel end
abstract type CPCRDFModel <: EoSModel end

struct CubicPlusChain{M} <: CubicPlusChainModel
    components::Vector{String}
    cubicmodel::M
    references::Vector{String}
end

include("chainmodel/chainmodel.jl")
include("cubicmodel/cubicmodel.jl")

function CubicPlusChain(components::Vector{String}; 
    idealmodel=BasicIdeal,
    cubicmodel = RK,
    alpha,
    mixing,
    rdf = ElliottRDF,
    translation=NoTranslation,
    userlocations=String[], 
    ideal_userlocations=String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    rdf_userlocations = String[],
    translation_userlocations = String[],
    verbose=false)

    #hack: we supplant activity with rdf
    activity = rdf
    activity_userlocations = rdf_userlocations

    cubic = cubicmodel(components;idealmodel,alpha,mixing,translation,activity,userlocations,ideal_userlocations,alpha_userlocations,mixing_userlocations,translation_userlocations,activity_userlocations,verbose)
    references = ["10.1021/acs.iecr.9b00435","10.1021/acs.iecr.9b00436","10.1021/acs.iecr.0c02483"]
    return CubicPlusChain(components,cubic,references)
end

idealmodel(model::CubicPlusChainModel) = idealmodel(model.cubicmodel)

function a_res(model::CubicPlusChainModel,V,T,z,_data = @f(data))
    cubicdata,m̄,β = _data
    return m̄*a_res(model.cubicmodel,V,T,z,cubicdata) + a_chain(model,V,T,z,_data)
end

function g_rdf end

function a_chain(model::CubicPlusChainModel,V,T,z,_data = @f(data))
    cubicdata,m̄,β = data
    RDF = model.cubicmodel.mixing.rdf
    return (m̄ - 1)*log(g_rdf(RDF,β))
end

function data(model::CubicPlusChainModel,V,T,z)
    m = model.cubicmodel.mixing.params.segment.values
    m̄ = dot(m,z)/∑z
    cubicdata = data(model.cubicmodel,V,T,z)
    ∑z,a,b,c = cubicdata
    β = ∑z*m̄*b/V
    return cubicdata,m̄,β
end

include("variants/CPC_RKE.jl")
include("variants/CPC_SRKE.jl")
include("variants/CPC_SRKE_bT.jl")


