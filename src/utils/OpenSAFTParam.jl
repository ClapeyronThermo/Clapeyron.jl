abstract type OpenSAFTParam end

struct SingleParam{T} <: OpenSAFTParam
    name::String
    values::Array{T,1}
    ismissingvalues::Array{Bool,1}
    components::Array{String,1}
    allcomponentsites::Array{Array{String,1},1}
    sourcecsvs::Array{String,1}
    sources::Array{String,1}
end

function SingleParam(x::SingleParam{T}) where T
    return SingleParam(x.name, deepcopy(x.values), deepcopy(x.ismissingvalues), x.components, x.allcomponentsites, x.sourcecsvs, deepcopy(x.sources))
end
function SingleParam(x::SingleParam{T}, v::Array{T,1}) where T
    return SingleParam(x.name, v, deepcopy(x.ismissingvalues), x.components, x.allcomponentsites, x.sourcecsvs, deepcopy(x.sources))
end

struct PairParam{T} <: OpenSAFTParam
    name::String
    values::Array{T,2}
    diagvalues::SubArray{T,1,Array{T,1},Tuple{Array{Int64,1}},false}
    ismissingvalues::Array{Bool,2}
    components::Array{String,1}
    allcomponentsites::Array{Array{String,1},1}
    sourcecsvs::Array{String,1}
    sources::Array{String,1}
    function PairParam(name::String, values::Array{T,2}, ismissingvalues::Array{Bool,2}, components::Array{String,1}, allcomponentsites::Array{Array{String,1},1}, sourcecsvs::Array{String,1}, sources::Array{String,1}) where T
        diagvalues = view(values, [1+(length(components)+1)*x for x in 0:length(components)-1])
        return new{T}(name, values, diagvalues, ismissingvalues, components, allcomponentsites, sourcecsvs, sources)
    end
end

function PairParam(x::PairParam{T}) where T
    return PairParam(x.name, deepcopy(x.values), deepcopy(x.ismissingvalues), x.components, x.allcomponentsites, x.sourcecsvs, deepcopy(x.sources))
end
function PairParam(x::PairParam{T}, v::Array{T,2}) where T
    return PairParam(x.name, v, deepcopy(x.ismissingvalues), x.components, x.allcomponentsites, x.sourcecsvs, deepcopy(x.sources))
end
function PairParam(x::SingleParam{T}) where T
    return PairParam(x.name, convertsingletopair(x.values), convert(Array{Bool},.!(convertsingletopair(convert(Array{Bool},.!(x.ismissingvalues))))), x.components, x.allcomponentsites, x.sourcecsvs, deepcopy(x.sources))
end
function PairParam(x::SingleParam{T}, v::Array{T,1}) where T
    return PairParam(x.name, convertsingletopair(v), convert(Array{Bool},.!(convertsingletopair(convert(Array{Bool},.!(x.ismissingvalues))))), x.components, x.allcomponentsites, x.sourcecsvs, deepcopy(x.sources))
end
function PairParam(x::SingleParam{T}, v::Array{T,2}) where T
    return PairParam(x.name, v, convert(Array{Bool},.!(convertsingletopair(convert(Array{Bool},.!(x.ismissingvalues))))), x.components, x.allcomponentsites, x.sourcecsvs, deepcopy(x.sources))
end

struct AssocParam{T} <: OpenSAFTParam
    name::String
    values::Array{Array{T,2},2}
    ismissingvalues::Array{Array{Bool,2},2}
    components::Array{String,1}
    allcomponentsites::Array{Array{String,1},1}
    sourcecsvs::Array{String,1}
    sources::Array{String,1}
end

function AssocParam(x::AssocParam{T}) where T
    return PairParam{T}(x.name, deepcopy(x.values), deepcopy(x.ismissingvalues), x.components, x.allcomponentsites, x.sourcecsvs, deepcopy(x.sources))
end
function AssocParam(x::AssocParam{T}, v::Array{Array{T,2},2}) where T
    return PairParam{T}(x.name, v, deepcopy(x.ismissingvalues), x.components, x.allcomponentsites, x.sourcecsvs, deepcopy(x.sources))
end

struct GCParam <: OpenSAFTParam
    components::Array{String,1}
    allcomponentgroups::Array{Array{String,1},1}
    allcomponentngroups::Array{Array{Int,1},1}
    flattenedgroups::Array{String,1}
    allcomponentnflattenedgroups::Array{Array{Int,1},1}
    sourcecsvs::Array{String,1}
end

struct SiteParam <: OpenSAFTParam
    components::Array{String,1}
    allcomponentsites::Array{Array{String,1},1}
    allcomponentnsites::Array{Array{Int,1},1}
    sourcecsvs::Array{String,1}
end
