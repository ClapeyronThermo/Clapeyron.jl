struct ElliottRDF <: CPCRDFModel
end

"""
    ElliottRDF <: CPCRDFModel
    ElliottRDF()

## Input parameters

None

## Description

Elliott Radial distance function.
```
    g_rdf(::ElliottRDF,β) = 1/(1 - 0.475*β)
```
"""
ElliottRDF

export ElliottRDF
function ElliottRDF(components::Array{String,1}; userlocations::Array{String,1}=String[], verbose=false)
    return ElliottRDF()
end

function ElliottRDF(; userlocations::Array{String,1}=String[], verbose=false)
    return ElliottRDF()
end

is_splittable(::ElliottRDF) = false

g_rdf(::ElliottRDF,β) = 1/(1 - 0.475*β)
