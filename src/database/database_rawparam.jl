#=This is just a tape.
component info stores (comp1,comp2,site1,site2). on singleparams, comp = comp1
on single and pair params, site1, site2, = ""
data is a vector of data found
sources is the sources for each point, same with csv
the strategy is to check each csv and produce RawParams, with ONLY nonmissing values
then we join the same params (joindata!) and finally, we "compile" the tapes (via compile_param)
we calculate the sites with the parsed raw params, as they have all the necessary information
For Clapeyron 0.4.0, this will also hold the group type, tapes with different group types cannot be merged.
=#
struct RawParam{T}
    name::String
    component_info::Vector{NTuple{4,String}}  # "Tape" for component data (component1,component2,site1,site2)
    data::Vector{T} # "Tape" of data
    sources::Vector{String} # "Tape" of parsed sources
    csv::Vector{String} # "Tape" of origin csv
    type::CSVType #the type of data
    grouptype::Symbol #in the future, this could hold an "Options" type,generated per CSV
end

Base.eltype(raw::RawParam) = Base.eltype(raw.data)
Base.length(raw::RawParam) = Base.length(raw.data)

function Base.show(io::IO,param::RawParam)
    print(io,typeof(param))
    print(io,param.component_info,",")
    print(io,param.data,",")
    print(io,param.type,")")
end

#join CSV types
#single and pair data can be merged onto pairdata
#other csv types cannot be merged.
# returns (type::CSVType, success::Bool)
function joindata!(old::CSVType,new::CSVType)
    if new == old
        return (new,true)
    elseif old ∈ (singledata,pairdata) && new ∈ (singledata,pairdata)
        return (pairdata,true)
    else
        return (old,false)
    end
end

#join two tapes
#tapes are destroyed here. in the sense that the "old" and "new" values are not valid anymore.
Base.@nospecialize

function joindata!(old::Vector,new::Vector)
    T1,T2 = eltype(old),eltype(new)
    if promote_type(T1,T2) == T1
        return append!(old,new)
    elseif promote_type(T1,T2) == T2
        return prepend!(new,old)
    else
        return vcat(string.(old),string.(new))
    end
end

function joindata!(old::RawParam,new::RawParam)
    tnew,type_sucess = joindata!(old.type,new.type)
    if old.grouptype !== new.grouptype
        if new.grouptype != :unkwown #for backwards compatibility
            error_different_grouptype(old,new)
        end
    end
    
    if !type_sucess
        error_clashing_headers(old,new)
    end
    component_info = append!(old.component_info,new.component_info)

    #Handle all the type variability of the data here
    data = joindata!(old.data,new.data)

    sources = append!(old.sources,new.sources)
    csv = append!(old.csv,new.csv)

    return RawParam(old.name,component_info,data,sources,csv,tnew,old.grouptype)
end

error_different_grouptype(old::RawParam,new::RawParam) = error_different_grouptype(old.grouptype,new.grouptype)

@noinline function error_different_grouptype(old::Symbol,new::Symbol)
    throw(error("""cannot join two databases with different group types:
    current group type: $(old)
    incoming group type: $(new)
    """))
end



@noinline function error_clashing_headers(old::RawParam,new::RawParam)
    told = Symbol(old.type)
    tnew = Symbol(new.type)
    header = error_color(old.name)
    err = """cannot join CSV header $header with incompatible data:
    current data type: $(told), with tables:
    - $(old.csv)
    incoming data type: $(tnew), with tables:
    - $(new.csv)
    """
    throw(error(err))
end

Base.@specialize

@noinline function error_clashing_headers(old::CSVType,new::CSVType,header)
    header = error_color(header)
    ("Header ", header, " appear ∈ both loaded assoc and non-assoc files.")
end
#=
compile_param takes a RawParam, and instantiates it into a ClapeyronParam.
it also builds empty params, if you pass a CSVType instead of a RawParam
=#

Base.@nospecialize

function compile_param(components,name,raw::RawParam,site_strings,options)
    if raw.type == singledata || raw.type == groupdata
        return compile_single(name,components,raw,options)
    elseif raw.type == pairdata
        return compile_pair(name,components,raw,options)
    elseif raw.type == assocdata
        return compile_assoc(name,components,raw,site_strings,options)
    end
    return nothing
end

function compile_param(components,name,raw::CSVType,site_strings,options)
    if raw == singledata
        return compile_single(name,components,raw,options)
    elseif raw == pairdata
        return compile_pair(name,components,raw,options)
    elseif raw == assocdata
        return compile_assoc(name,components,raw,site_strings,options)
    end
    return nothing
end

function compile_single(name,components,raw::RawParam,options)
    EMPTY_STR = ""
    l = length(components)
    L = eltype(raw)
    if L <: Number
        values = zeros(L,l)
    else
        values = fill("",l)
    end
    ismissingvals = ones(Bool,l)
    sources = fill(EMPTY_STR,l)
    sources_csv = fill(EMPTY_STR,l)
    for (k,v,ss,sc) ∈ zip(raw.component_info,raw.data,raw.sources,raw.csv)
        i = findfirst(==(k[1]),components)
        values[i] = v
        ismissingvals[i] = false
        sources[i] = ss
        sources_csv[i] = sc
    end
    sources = unique!(sources)
    sources_csv = unique!(sources_csv)
    filter!(!isequal(EMPTY_STR),sources)
    filter!(!isequal(EMPTY_STR),sources_csv)
    return SingleParameter(name,components,values,ismissingvals,sources_csv,sources)
end

function compile_single(name,components,type::CSVType,options)
    if name ∈ options.ignore_missing_singleparams
        return SingleParam(name,components)
    else
        error("cannot found any values for ", error_color(name), ".")
    end
end

function compile_pair(name,components,raw::RawParam,options)
    EMPTY_STR = ""
    symmetric = name ∉ options.asymmetricparams
    l = length(components)
    L = eltype(raw)
    if L <: Number
        values = zeros(L,(l,l))
    else
        values = fill("",(l,l))
    end
    ismissingvals = ones(Bool,(l,l))
    sources = fill(EMPTY_STR,(l,l))
    sources_csv = fill(EMPTY_STR,(l,l))
    for (k,v,ss,sc) ∈ zip(raw.component_info,raw.data,raw.sources,raw.csv)
        c1,c2,_,_ = k
        i = findfirst(==(c1),components)
        #if the second component is null, it comes from a single param, then i = (i,i)
        j = k[2] == "" ? i : findfirst(==(c2),components)
        values[i,j] = v
        ismissingvals[i,j] = false
        sources[i,j] = ss
        sources_csv[i,j] = sc
        if symmetric && i ≠ j
            ismissingvals[j,i] = false
            values[j,i] = v
            sources[j,i] = ss
            sources_csv[j,i] = sc
        end
    end
    sources = unique!(vec(sources))
    sources_csv = unique!(vec(sources_csv))
    filter!(!isequal(EMPTY_STR),sources)
    filter!(!isequal(EMPTY_STR),sources_csv)
    return PairParameter(name,components,values,ismissingvals,sources_csv,sources)
end

function compile_pair(name,components,type::CSVType,options)
    return PairParam(name,components)
end

function compile_assoc(name,components,raw::RawParam,site_strings,options)
    EMPTY_STR = ""
    _ijab = standarize_comp_info(raw.component_info,components,site_strings)
    unique_sitepairs = unique(raw.component_info)
    l = length(unique_sitepairs)
    unique_dict = Dict{NTuple{4,String},Int}(unique_sitepairs[i] => i for i in 1:l)
    sources_csv = fill(EMPTY_STR,l)
    sources = fill(EMPTY_STR,l)
    ijab = similar(_ijab,l)
    inner_values = similar(raw.data,l)
    for (j,k) ∈ enumerate(raw.component_info)
        i = unique_dict[k]
        inner_values[i] = raw.data[j]
        ijab[i] = _ijab[j]
        sources[i] = raw.sources[j]
        sources_csv[i] = raw.csv[j]
    end
    
    idxs = sortperm(ijab) #CompressedAssoc4DMatrix requires lexicographically sorted component-site idxs
    ijab = ijab[idxs]
    inner_values = inner_values[idxs]
    sources = sources[idxs]
    sources_csv = sources_csv[idxs]
    values = Compressed4DMatrix(inner_values,ijab)
    unique!(sources)
    unique!(sources_csv)
    filter!(!isequal(EMPTY_STR),sources)
    filter!(!isequal(EMPTY_STR),sources_csv)
    param = AssocParam(name,components,values,site_strings,sources_csv,sources)
    return param
end

function compile_assoc(name,components,raw::CSVType,site_strings,options)
    values = Compressed4DMatrix{Float64}()
    return AssocParam(name,components,values,site_strings,String[],String[])
end

#Sort site tape, so that components are sorted by the input.
function standarize_comp_info(component_info,components,site_strings)
    ijab = Vector{Tuple{Int,Int,Int,Int}}(undef,length(component_info))
    l = length(components)
    for (i,val) ∈ pairs(component_info)
        c1,c2,s1,s2 = val
        idx1 = findfirst(isequal(c1), components)
        idx2 = findfirst(isequal(c2), components)
        idx21 = findfirst(isequal(s1), site_strings[idx1])
        idx22 = findfirst(isequal(s2), site_strings[idx2])
        if idx1 > idx2
            newval = (c2,c1,s2,s1)
            ijab[i] = (idx2,idx1,idx22,idx21)
        elseif idx1 < idx2
            newval = val
            ijab[i] = (idx1,idx2,idx21,idx22)
        else
            if idx21 > idx22
                ijab[i] = (idx1,idx2,idx22,idx21)
                newval = (c1,c1,s2,s1)
            else
                ijab[i] = (idx1,idx2,idx21,idx22)
                newval = val
            end
        end
        component_info[i] = newval
    end
    return ijab
end
#=check valid params
For single params, it checks that there aren't missing values (can be overrided)
For pair params, it checks that there aren't missing values ∈ the diagonal.  one exception is when all
values are zero, ∈ this case is ommited (can also be overrided)
=#
function is_valid_param(param::SingleParameter,options)
    missingvals = param.ismissingvalues
    if param.name ∉ options.ignore_missing_singleparams && any(missingvals)
        vals = [ifelse(missingvals[i],missing,param.values[i]) for i ∈ 1:length(missingvals)]
        error("Missing values exist ∈ single parameter ", error_color(param.name), ": ", vals, ".")
    end
    return nothing
end

function is_valid_param(param::PairParameter,options)
    diag = diagvalues(param.ismissingvalues)
    if param.name ∉ options.ignore_missing_singleparams && !all(diag) && any(diag)
        vals = [ifelse(diag[i],missing,param.values[i]) for i ∈ 1:length(diag)]
        error("Partial missing values exist ∈ diagonal of pair parameter ",error_color(param.name), ": ", vals, ".")
    end
    return nothing
end

function is_valid_param(param::AssocParam,options)
    return nothing
end

Base.@specialize