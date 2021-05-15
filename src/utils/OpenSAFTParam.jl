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

function Base.show(io::IO,param::SingleParam)
    print(io,"SingleParam(")
    print(io,"\"",param.name,"\"",")[")
    for (name,val,miss,i) in zip(param.components,param.values,param.ismissingvalues,1:length(param.values))
        i != 1 && print(io,",")
        if miss == false
            print(io,name,"=",val)
        else
            print(io,name,"=","-")
        end
    end
    print(io,"]")
end

function Base.show(io::IO,::MIME"text/plain",param::SingleParam)
    len = length(param.values)
    println(io,"SingleParam(", "\"",param.name,"\")"," with ",len," component",ifelse(len==1,":","s:"))
    i = 0
    for (name,val,miss) in zip(param.components,param.values,param.ismissingvalues)
        i += 1
        if i > 1
            println(io)
        end
        if miss == false
            print(io," ",name," = ",val)
        else
            print(io," ",name," = ","-")
        end
    end
end



function SingleParam(x::SingleParam{T}) where T
    return SingleParam(x.name, deepcopy(x.values), deepcopy(x.ismissingvalues), x.components, x.allcomponentsites, deepcopy(x.sourcecsvs), deepcopy(x.sources))
end

#a barebones constructor, in case we dont build from csv
function SingleParam(
    name::String,
    components::Vector{String},
    values::Vector{T},
    ismissingvalues::Vector{Bool} = [ismissing for i = 1:length(value)],
    allcomponentsites = Array{Array{String,1},1}(undef,0),
    sourcecsvs = String[],
    sources = String[]
    ) where T<:Union{ <: Real,Missing}
    
    _values,_ismissingvalues = nondefaultmissing(values; defaultvalue=nothing)
    TT = eltype(_values)
    return  SingleParam{TT}(names, _values, _ismissingvalues, components, allcomponentsites, sourcecsvs, sources)
end


function SingleParam(x::SingleParam, v::Array{T,1}) where T
    return SingleParam(x.name, v, deepcopy(x.ismissingvalues), x.components, x.allcomponentsites, deepcopy(x.sourcecsvs), deepcopy(x.sources),)
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
    return PairParam(x.name, deepcopy(x.values), deepcopy(x.ismissingvalues), x.components, x.allcomponentsites, deepcopy(x.sourcecsvs), deepcopy(x.sources))
end
function PairParam(x::SingleParam{T}) where T
    return PairParam(x.name, convertsingletopair(x.values), convert(Array{Bool},.!(convertsingletopair(convert(Array{Bool},.!(x.ismissingvalues))))), x.components, x.allcomponentsites, deepcopy(x.sourcecsvs), deepcopy(x.sources))
end

function PairParam(x::PairParam, v::Array{T,2}) where T
    return PairParam(x.name, v, deepcopy(x.ismissingvalues), x.components, x.allcomponentsites, deepcopy(x.sourcecsvs), deepcopy(x.sources))
end
function PairParam(x::SingleParam, v::Array{T,1}) where T
    return PairParam(x.name, convertsingletopair(v), convert(Array{Bool},.!(convertsingletopair(convert(Array{Bool},.!(x.ismissingvalues))))), x.components, x.allcomponentsites, deepcopy(x.sourcecsvs), deepcopy(x.sources))
end
function PairParam(x::SingleParam, v::Array{T,2}) where T
    return PairParam(x.name, v, convert(Array{Bool},.!(convertsingletopair(convert(Array{Bool},.!(x.ismissingvalues))))), x.components, x.allcomponentsites, deepcopy(x.sourcecsvs), deepcopy(x.sources))
end

#=
function Base.show(io::IO,::MIME"text/plain",param::PairParam)
    len = length(param.values)
    println(io,"PairParam(", "\"",param.name,"\")"," with ",len," component",ifelse(len==1,":","s:"))
    i = 0
    for (name,val,miss) in zip(param.components,param.values,param.ismissingvalues)
        i += 1
        if i > 1
            println(io)
        end
        if miss == false
            print(io," ",name," = ",val)
        else
            print(io," ",name," = ","-")
        end
    end
end

function Base.show(io::IO,param::PairParam)
    print(io,"PairParam(")
    print(io,"\"",param.name,"\"",")[")
    for (name,val,miss,i) in zip(param.components,param.values,param.ismissingvalues,1:length(param.values))
        i != 1 && print(io,",")
        if miss == false
            print(io,name,"=",val)
        else
            print(io,name,"=","-")
        end
    end
    print(io,"]")
end
=#
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
    return PairParam{T}(x.name, deepcopy(x.values), deepcopy(x.ismissingvalues), x.components, x.allcomponentsites, deepcopy(x.sourcecsvs), deepcopy(x.sources))
end

function AssocParam{T}(x::AssocParam, v::Array{Array{T,2},2}) where T
    return AssocParam{T}(x.name, v, deepcopy(x.ismissingvalues), x.components, x.allcomponentsites, deepcopy(x.sourcecsvs), deepcopy(x.sources))
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

function SiteParam(pairs::Dict{String,SingleParam{Int}})
    arbitraryparam = first(values(pairs))
    components = arbitraryparam.components
    allcomponentsites = arbitraryparam.allcomponentsites
    sourcecsvs = unique([([x.sourcecsvs for x in values(pairs)]...)...])
    allcomponentnsites = [[pairs[allcomponentsites[i][j]].values[i] for j ∈ 1:length(allcomponentsites[i])] for i ∈ 1:length(components)]  # or groupsites
    return SiteParam(components, allcomponentsites, allcomponentnsites, sourcecsvs)
end

paramvals(param::OpenSAFTParam) = param.values
paramvals(x) = x
