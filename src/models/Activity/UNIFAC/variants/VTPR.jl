struct VTPRUNIFACParam <: EoSParam
    A::PairParam{Float64}
    B::PairParam{Float64}
    C::PairParam{Float64}
    Q::SingleParam{Float64}
end

abstract type VTPRUNIFACModel <: UNIFACModel end

struct VTPRUNIFAC{c<:EoSModel} <: VTPRUNIFACModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    groups::GroupParam
    params::VTPRUNIFACParam
    puremodel::Vector{c}
    absolutetolerance::Float64
    references::Array{String,1}
end

@registermodel VTPRUNIFAC
export VTPRUNIFAC

function VTPRUNIFAC(components; puremodel=PR,
    userlocations=String[], 
     verbose=false)
    groups = GroupParam(components, ["Activity/UNIFAC/VTPR/VTPR_groups.csv"]; verbose=verbose)

    params = getparams(groups, ["Activity/UNIFAC/VTPR/VTPR_like.csv", "Activity/UNIFAC/VTPR/VTPR_unlike.csv"]; userlocations=userlocations,  asymmetricparams=["A","B","C"], ignore_missing_singleparams=["A","B","C"], verbose=verbose)
    A  = params["A"]
    B  = params["B"]
    C  = params["C"]
    Q  = params["Q"]
    icomponents = 1:length(components)
    
    init_puremodel = [puremodel([groups.components[i]]) for i in icomponents]
    packagedparams = VTPRUNIFACParam(A,B,C,Q)
    references = String[]
    model = VTPRUNIFAC(components,icomponents,groups,packagedparams,init_puremodel,1e-12,references)
    return model
end

function activity_coefficient(model::VTPRUNIFACModel,V,T,z)
    return exp.(@f(lnγ_res))
end