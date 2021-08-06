abstract type ClapeyronParam end

struct SingleParam{T} <: ClapeyronParam
    name::String
    values::Array{T,1}
    ismissingvalues::Array{Bool,1}
    components::Array{String,1}
    allcomponentsites::Array{Array{String,1},1}
    sourcecsvs::Array{String,1}
    sources::Array{String,1}
end

function Base.show(io::IO, param::SingleParam)
    print(io, typeof(param), "(\"", param.name, "\")[")
    for component in param.components
        component != first(param.components) && print(io, ",")
        print(io, "\"", component, "\"")
    end
    print(io, "]")
end

function Base.show(io::IO, ::MIME"text/plain", param::SingleParam)
    len = length(param.values)
    print(io, typeof(param), "(\"", param.name)
    println(io, "\") with ", len, " component", ifelse(len==1, ":", "s:"))
    i = 0
    for (name, val, miss) in zip(param.components, param.values, param.ismissingvalues)
        i += 1
        if i > 1
            println(io)
        end
        if miss == false
            if typeof(val) <: AbstractString
                print(io, " \"", name, "\" => \"", val, "\"")
            else
                print(io, " \"", name, "\" => ", val)
            end
        else
            print(io, " \"", name, " => -")
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

struct PairParam{T} <: ClapeyronParam
    name::String
    values::Array{T,2}
    diagvalues::SubArray{T, 1, Vector{T}, Tuple{StepRange{Int64, Int64}}, true}
    ismissingvalues::Array{Bool,2}
    components::Array{String,1}
    allcomponentsites::Array{Array{String,1},1}
    sourcecsvs::Array{String,1}
    sources::Array{String,1}
end

function PairParam(name::String, values::Array{T,2}, ismissingvalues::Array{Bool,2}, components::Array{String,1}, allcomponentsites::Array{Array{String,1},1}, sourcecsvs::Array{String,1}, sources::Array{String,1}) where T
    diagvalues = view(values, diagind(values))
    return PairParam{T}(name, values, diagvalues, ismissingvalues, components, allcomponentsites, sourcecsvs, sources)
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

function Base.show(io::IO,mime::MIME"text/plain",param::PairParam{T}) where T
    print(io,"PairParam{",string(T),"}")
    show(io,param.components)
    println(io,") with values:")
    show(io,mime,param.values)
end

function Base.show(io::IO,param::PairParam)
    print(io,"PairParam(",param.name)
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

struct AssocParam{T} <: ClapeyronParam
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

function Base.show(io::IO,mime::MIME"text/plain",param::AssocParam{T}) where T
    print(io,"AssocParam{",string(T),"}")
    show(io,param.components)
    println(io,") with values:")
    show(io,mime,param.values)
end

function Base.show(io::IO,param::AssocParam)
    print(io,"AssocParam(",param.name)
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

struct GroupParam <: ClapeyronParam
    components::Array{String,1}
    groups::Array{Array{String,1},1}
    n_groups::Array{Array{Int,1},1}
    i_groups::Array{Array{Int,1},1}
    flattenedgroups::Array{String,1}
    n_flattenedgroups::Array{Array{Int,1},1}
    i_flattenedgroups::UnitRange{Int}
    sourcecsvs::Array{String,1}
end

function Base.show(io::IO, mime::MIME"text/plain", param::GroupParam)
    print(io,"GroupParam ")
    len = length(param.components)
    println(io,"with ", len, " component", ifelse(len==1, ":", "s:"))
    
    for i in 1:length(param.components)
        
        print(io, " \"", param.components[i], "\": ")
        firstloop = true
        for j in 1:length(param.n_groups[i])
            firstloop == false && print(io, ", ")
            print(io, "\"", param.groups[i][j], "\" => ", param.n_groups[i][j])
            firstloop = false
        end
        i != length(param.components) && println(io)
    end 
end

function Base.show(io::IO, param::GroupParam)
    print(io,"GroupParam[")
    len = length(param.components)
    
    for i in 1:length(param.components)
        
        print(io, "\"", param.components[i], "\" => [")
        firstloop = true
        for j in 1:length(param.n_groups[i])
            firstloop == false && print(io, ", ")
            print(io, "\"", param.groups[i][j], "\" => ", param.n_groups[i][j])
            firstloop = false
        end
        print(io,']')
        i != length(param.components) && print(io,", ")
    end
    print(io,"]")
end


struct SiteParam <: ClapeyronParam
    components::Array{String,1}
    sites::Array{Array{String,1},1}
    n_sites::Array{Array{Int,1},1}
    length_sites::Array{Int,1}
    i_sites::Array{UnitRange{Int},1}
    sourcecsvs::Array{String,1}
end

function Base.show(io::IO, mime::MIME"text/plain", param::SiteParam)
    print(io,"SiteParam ")
    len = length(param.components)
    println(io,"with ", len, " site", ifelse(len==1, ":", "s:"))
    
    for i in 1:length(param.components)
        
        print(io, " \"", param.components[i], "\": ")
        firstloop = true
        if length(param.n_sites[i]) == 0
            print(io,"(no sites)")
        end
        for j in 1:length(param.n_sites[i])
            firstloop == false && print(io, ", ")
            print(io, "\"", param.sites[i][j], "\" => ", param.n_sites[i][j])
            firstloop = false
        end
        i != length(param.components) && println(io)
    end 
end

function Base.show(io::IO, param::SiteParam)
    print(io,"SiteParam[")
    len = length(param.components)
    
    for i in 1:length(param.components)
        
        print(io, "\"", param.components[i], "\" => [")
        firstloop = true
    
        for j in 1:length(param.n_sites[i])
            firstloop == false && print(io, ", ")
            print(io, "\"", param.sites[i][j], "\" => ", param.n_sites[i][j])
            firstloop = false
        end
        print(io,']')
        i != length(param.components) && print(io,", ")
    end
    print(io,"]")
end


function SiteParam(pairs::Dict{String,SingleParam{Int}})
    arbitraryparam = first(values(pairs))
    components = arbitraryparam.components
    sites = arbitraryparam.allcomponentsites
    sourcecsvs = unique([([x.sourcecsvs for x in values(pairs)]...)...])
    n_sites = [[pairs[sites[i][j]].values[i] for j ∈ 1:length(sites[i])] for i ∈ 1:length(components)]  # or groupsites
    length_sites = [length(componentsites) for componentsites ∈ sites]
    i_sites = [1:length_sites[i] for i ∈ 1:length(components)]

    return SiteParam(components,
    sites,
    n_sites,
    length_sites,
    i_sites,
    sourcecsvs)
end

#empty SiteParam
function SiteParam(components::Vector{String})
    n = length(components)
    return SiteParam(components,
    [String[] for _ ∈ 1:n],
    [Int[] for _ ∈ 1:n], 
    zeros(Int,n),
    [1:0 for _ ∈ 1:n],
    String[])
end

paramvals(param::ClapeyronParam) = param.values
paramvals(x) = x


function split_model(param::SingleParam{T},components = param.components) where T
    return [SingleParam{T}(
    param.name,
    [param.values[i]],
    [param.ismissingvalues[i]],
    [components[i]],
    [param.allcomponentsites[i]],
    deepcopy(param.sourcecsvs),
    deepcopy(param.sources)
    ) for i in 1:length(components)]
end

#this conversion is lossy, as interaction between two or more components are lost.
function split_model(param::PairParam{T},components = param.components) where T
    function generator(i)     
        _value = zeros(T,1,1)
        _value[1,1] = param.values[i,i]
        _diagvalue = view(_value,1:1:1)
        _ismissingvalues = zeros(Bool,1,1)
        _ismissingvalues[1,1] = param.ismissingvalues[i,i]
        return PairParam{T}(
                param.name,
                _value,
                _diagvalue,
                _ismissingvalues,
                [components[i]],
                [param.allcomponentsites[i]],
                deepcopy(param.sourcecsvs),
                deepcopy(param.sources)
                )
    end 
    return [generator(i) for i in 1:length(components)]
end

function split_model(param::Vector,components)
    checkbounds(length(param),length(components))
    return [[xi] for xi in param]
    
end

#this conversion is lossy, as interaction between two or more components are lost.
function split_model(param::AssocParam{T},components = param.components) where T
    function generator(i)     
        _value = Matrix{Matrix{T}}(undef,1,1)
        _value[1,1] = param.values[i,i]
        _ismissingvalues = Matrix{Matrix{Bool}}(undef,1,1)
        _ismissingvalues[1,1] = param.ismissingvalues[i,i]
        return AssocParam{T}(
                param.name,
                _value,
                _ismissingvalues,
                [components[i]],
                [param.allcomponentsites[i]],
                deepcopy(param.sourcecsvs),
                deepcopy(param.sources)
                )
    end 
    return [generator(i) for i in 1:length(components)]
end

export SingleParam, SiteParam, PairParam, AssocParam, GroupParam
