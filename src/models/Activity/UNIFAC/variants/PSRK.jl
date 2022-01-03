
abstract type PSRKUNIFACModel <: UNIFACModel end

struct PSRKUNIFAC{c<:EoSModel} <: PSRKUNIFACModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    groups::GroupParam
    params::UNIFACParam
    puremodel::Vector{c}
    absolutetolerance::Float64
    references::Array{String,1}
end

@registermodel PSRKUNIFAC
export PSRKUNIFAC

function PSRKUNIFAC(components; puremodel=PR,
    userlocations=String[], 
     verbose=false)
    groups = GroupParam(components, ["Activity/UNIFAC/PSRK/PSRK_groups.csv"]; verbose=verbose)

    params = getparams(groups, ["Activity/UNIFAC/PSRK/PSRK_like.csv", "Activity/UNIFAC/PSRK/PSRK_unlike.csv"]; userlocations=userlocations,  asymmetricparams=["A","B","C"], ignore_missing_singleparams=["A","B","C"], verbose=verbose)
    A  = params["A"]
    B  = params["B"]
    C  = params["C"]
    R  = params["R"]
    Q  = params["Q"]
    icomponents = 1:length(components)
    
    init_puremodel = [puremodel([groups.components[i]]) for i in icomponents]
    packagedparams = UNIFACParam(A,B,C,R,Q)
    references = String[]
    model = PSRKUNIFAC(components,icomponents,groups,packagedparams,init_puremodel,1e-12,references)
    return model
end
