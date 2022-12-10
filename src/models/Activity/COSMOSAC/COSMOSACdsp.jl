struct COSMOSACdspParam <: EoSParam
    Pnhb::SingleParam{Vector{Float64}}
    POH::SingleParam{Vector{Float64}}
    POT::SingleParam{Vector{Float64}}
    epsilon::SingleParam{Float64}
    V::SingleParam{Float64}
    A::SingleParam{Float64}
    water::SingleParam{Int64}
    COOH::SingleParam{Int64}
    hb_acc::SingleParam{Int64}
    hb_don::SingleParam{Int64}
end

abstract type COSMOSACdspModel <: COSMOSAC10Model end

struct COSMOSACdsp{c<:EoSModel} <: COSMOSACdspModel
    components::Array{String,1}
    params::COSMOSACdspParam
    puremodel::EoSVectorParam{c}
    absolutetolerance::Float64
    references::Array{String,1}
end

@registermodel COSMOSACdsp
export COSMOSACdsp

function COSMOSACdsp(components::Vector{String};
    puremodel = PR,
    userlocations = String[],
    pure_userlocations = String[],
    verbose=false, kwargs...)

    params = getparams(components, ["Activity/COSMOSAC/COSMOSAC10_like.csv","Activity/COSMOSAC/COSMOSACdsp_like.csv"]; userlocations=userlocations, verbose=verbose)
    Pnhb  = COSMO_parse_Pi(params["Pnhb"])
    POH  = COSMO_parse_Pi(params["POH"])
    POT  = COSMO_parse_Pi(params["POT"])
    A  = params["A"]
    V  = params["V"]
    epsilon  = params["epsilon"]
    water = params["water"]
    COOH = params["COOH"]
    hb_acc = params["hb_acc"]
    hb_don = params["hb_don"]
    
    _puremodel = init_puremodel(puremodel,components,pure_userlocations,verbose)
    packagedparams = COSMOSACdspParam(Pnhb,POH,POT,epsilon,V,A,water,COOH,hb_acc,hb_don)
    references = String[]
    model = COSMOSACdsp(components,packagedparams,_puremodel,1e-12,references)
    return model
end

function activity_coefficient(model::COSMOSACdspModel,V,T,z)
    return exp.(@f(lnγ_comb) .+@f(lnγ_res).+@f(lnγ_dsp))
end

function lnγ_dsp(model::COSMOSACdspModel,V,T,z)
    x = z./sum(z)
    if model.params.water.values[1]==1 && model.params.hb_acc.values[2]==1 && model.params.hb_don.values[2]==0
        w = -0.27027
    elseif model.params.water.values[2]==1 && model.params.hb_acc.values[1]==1 && model.params.hb_don.values[1]==0
        w = -0.27027
    elseif model.params.COOH.values[1]==1 && (model.params.hb_acc.values[2]==1 | model.params.hb_don.values[2]==1)
        w = -0.27027
    elseif model.params.COOH.values[1]==1 && model.params.water.values[2]==1
        w = -0.27027
    elseif model.params.COOH.values[2]==1 && model.params.water.values[1]==1
        w = -0.27027
    else
        w = 0.27027
    end

    ϵ = model.params.epsilon.values

    A = w*(0.5*(ϵ[1]+ϵ[2])-√(ϵ[1]*ϵ[2]))
    return A*(1 .-x).^2
end
