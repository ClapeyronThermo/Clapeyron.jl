
abstract type HANNAModel <: ActivityModel end

struct HANNA{c<:EoSModel} <: HANNAModel
    components::Array{String,1}
    params::WilsonParam
    puremodel::EoSVectorParam{c}
    references::Array{String,1}
end

export HANNA

"""
    HANNA <: ActivityModel
    HANNA(components;
    puremodel = BasicIdeal(),
    userlocations = String[],
    pure_userlocations = String[],
    verbose = false)

## Input parameters

"""
HANNA

function HANNA(components;
        puremodel = BasicIdeal(),
        userlocations = String[],
        pure_userlocations = String[],
        verbose = false)
    
    return HANNA()
end

function excess_gibbs_free_energy(model::HANNAModel,p,T,z)
    
    return 
end
