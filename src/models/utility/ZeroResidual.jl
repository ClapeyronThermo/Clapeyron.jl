struct ZeroResidual{I} <: EoSModel
    components::Vector{String}
    idealmodel::I
end

Base.length(model::ZeroResidual) = length(model.components)
Base.show(io::IO,model::ZeroResidual) = eosshow(io,model)
Base.show(io::IO, mime::MIME"text/plain", model::ZeroResidual) = eosshow(io,mime,model)

"""
    ZeroResidual <: EoSModel
    ZeroResidual(components; 
    idealmodel=BasicIdeal,
    ideal_userlocations=String[],
    verbose=false)
## Input parameters
None
## Description
Zero residual model.
```
    a_res(model,V,T,z) = 0
```
"""
function ZeroResidual(components; 
    idealmodel=BasicIdeal,
    ideal_userlocations=String[],
    verbose=false,
    )
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    return ZeroResidual(components,init_idealmodel)
end

a_res(model::ZeroResidual,V,T,z) = 0.0
export ZeroResidual