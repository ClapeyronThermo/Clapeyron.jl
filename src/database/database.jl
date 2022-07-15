@enum CSVType singledata pairdata assocdata groupdata

#=This is just a tape.
component info stores (comp1,comp2,site1,site2). on singleparams, comp = comp1
on single and pair params, site1, site2, = ""
data is a vector of data found
sources is the sources for each point, same with csv

the strategy is to check each csv and produce RawParams, with ONLY nonmissing values
then we join the same params (joindata!) and then we "compile" the tapes (via compile_param)

we calculate the sites with the parsed raw params, as they have all the necessary information
For Clapeyron 0.4.0, this will also hold the group type, tapes with different group types cannot be merged.
=#
struct RawParam{T}
    name::String
    component_info::Vector{NTuple{4,String}}
    data::Vector{T}
    sources::Vector{String}
    csv::Vector{String}
    type::CSVType
    grouptype::String #in the future, this could hold an "Options" type,generated per CSV
end

Base.eltype(raw::RawParam) = Base.eltype(raw.data)

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
    elseif old in (singledata,pairdata) && new in (singledata,pairdata)
        return (pairdata,true)
    else
        return (old,false)
    end
end

#join two tapes
#tapes are destroyed here.
Base.@nospecialize

function joindata!(old::Vector,new::Vector)
    T1,T2 = eltype(old),eltype(new)
    if promote_type(T1,T2) == T1
        return prepend!(old,new)
    elseif promote_type(T1,T2) == T2
        return append!(new,old)
    else
        return vcat(string.(new),string.(old))
    end
end

function joindata!(old::RawParam,new::RawParam)
    component_info = prepend!(old.component_info,new.component_info)
    
    #Handle all the type variability of the data here
    data = joindata!(old.data,new.data)

    sources = prepend!(old.sources,new.sources)
    csv = prepend!(old.csv,new.csv)
    tnew,type_sucess = joindata!(old.type,new.type)
    if old.grouptype !== new.grouptype
        error_different_grouptype(old,new)
    end
    if !type_sucess
        error_clashing_headers(old,new)
    end
    return RawParam(old.name,component_info,data,sources,csv,tnew,old.grouptype)
end

Base.@specialize

@noinline function error_different_grouptype(old,new)
    throw(error("cannot join two databases with different group types! (old = $(old.grouptype), new = $(new.grouptype)"))
end

@noinline function error_clashing_headers(old::RawParam,new::RawParam)
    told = Symbol(old.type)
    tnew = Symbol(new.type)
    c = Base.text_colors
    reset = c[:normal]
    red = c[:bold] * c[:red]
    header = red * old.name * reset
    err = """cannot join CSV header $k with incompatible data:
    current data type: $(told), with tables:
    - $(old.csv)
    incoming data type: $(tnew), with tables:
    - $(new.csv)    
    """
    throw(error(err))
end

@noinline function error_clashing_headers(old::CSVType,new::CSVType,header)
    c = Base.text_colors
    reset = c[:normal]
    red = c[:bold] * c[:red]
    header = red * old.name * reset
    ("Header ", header, " appear in both loaded assoc and non-assoc files.")
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
    for (k,v,ss,sc) in zip(raw.component_info,raw.data,raw.sources,raw.csv)
        i = findfirst(==(k[1]),components)
        values[i] = v
        ismissingvals[i] = false
        sources[i] = ss
        sources_csv[i] = sc
    end
    sources = unique!(vec(sources))
    sources_csv = unique!(vec(sources_csv))
    return SingleParameter(name,components,values,ismissingvals,sources_csv,sources)
end

function compile_single(name,components,type::CSVType,options)
    l = length(components)
    values = zeros(Float64,l)
    ismissingvals = fill(true,l)
    return SingleParameter(name,components,values,ismissingvals,String[],String[])
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
    for (k,v,ss,sc) in zip(raw.component_info,raw.data,raw.sources,raw.csv)
        c1,c2,_,_ = k
        i = findfirst(==(c1),components)
        j = k[2] == "" ? i : findfirst(==(c2),components)
        values[i,j] = v
        ismissingvals[i,j] = false
        sources[i,j] = ss
        sources_csv[i,j] = sc
        if symmetric 
            values[j,i] = v
            ismissingvals[j,i] = false
            sources[j,i] = ss
            sources_csv[j,i] = sc
        end  
    end
    sources = unique!(vec(sources))
    sources_csv = unique!(vec(sources_csv))
    diagvals = diagvalues(values)
    return PairParameter(name,components,values,diagvals,ismissingvals,sources_csv,sources)
end

function compile_pair(name,components,type::CSVType,options)
    l = length(components)
    values = zeros(Float64,(l,l))
    ismissingvals = fill(true,(l,l))
    diagvals = diagvalues(values)
    return PairParameter(name,components,values,diagvals,ismissingvals,String[],String[])
end

function compile_assoc(name,components,raw::RawParam,site_strings,options)
    EMPTY_STR = ""
    vals = raw.data
    c_12,s_12 = standarize_comp_info(raw.component_info,components,site_strings)
    comp_info = raw.component_info
    idxs = sort!([last(findall(==(u),comp_info)) for u in unique(comp_info)])
    s_12 = s_12[idxs]
    c_12 = c_12[idxs]
    sources = raw.sources[idxs]
    csvs = raw.csv[idxs]
    ij = maximum(maximum(i) for i in c_12)
    ab = maximum(maximum(i) for i in s_12)
    size_ij = (ij,ij)
    size_ab =  (ab,ab)
    values = Compressed4DMatrix(vals[idxs],c_12,s_12,size_ij,size_ab)
    param = AssocParam(name,components,values,site_strings,csvs,sources)
    return param
end

function compile_assoc(name,components,raw::CSVType,site_strings,options)
    vals = Float64[]
    c_12 = Vector{Tuple{Int,Int}}(undef,0)
    s_12 = Vector{Tuple{Int,Int}}(undef,0)
    values = Compressed4DMatrix(vals,c_12,s_12,(0,0),(0,0))
    return AssocParam(name,components,values,site_strings,String[],String[])
end

#Sort site tape, so that components are sorted by the input.
function standarize_comp_info(component_info,components,site_strings)
    c_12 = Vector{Tuple{Int,Int}}(undef,length(component_info))
    s_12 = Vector{Tuple{Int,Int}}(undef,length(component_info))
    for (i,val) in pairs(component_info)
        c1,c2,s1,s2 = val
        idx1 = findfirst(isequal(c1), components)
        idx2 = findfirst(isequal(c2), components)
        idx21 = findfirst(isequal(s1), site_strings[idx1])
        idx22 = findfirst(isequal(s2), site_strings[idx2])
        if idx1 > idx2
            newval = (c2,c1,s2,s1)
            c_12[i] = (idx2,idx1)
            s_12[i] = (idx22,idx21)
        else
            newval = val
            c_12[i] = (idx1,idx2)
            s_12[i] = (idx21,idx22)
        end
        component_info[i] = newval
    end
    return c_12,s_12
end
#=check valid params
For single params, it checks that there aren't missing values (can be overrided)
For pair params, it checks that there aren't missing values in the diagonal.  one exception is when all
values are zero, in this case is ommited (can also be overrided)
=#
function is_valid_param(param::SingleParameter,options)
    missingvals = param.ismissingvalues
    if param.name ∉ options.ignore_missing_singleparams && any(missingvals)
        vals = [ifelse(missingvals[i],missing,param.values[i]) for i in 1:length(missingvals)]
        c = Base.text_colors
        reset = c[:normal]
        red = c[:bold] * c[:red]
        error("Missing values exist in single parameter ", red, param.name,reset, ": ", vals, ".")
    end
    return nothing
end

function is_valid_param(param::PairParameter,options)
    
    diag = diagvalues(param.ismissingvalues)
    if param.name ∉ options.ignore_missing_singleparams && !all(diag) && any(diag)
        c = Base.text_colors
        reset = c[:normal]
        red = c[:bold] * c[:red]
        vals = [ifelse(diag[i],missing,param.values[i]) for i in 1:length(diag)]
        error("Partial missing values exist in diagonal of pair parameter ",red,param.name,reset, ": ", vals, ".")
    end
    return nothing
end

function is_valid_param(param::AssocParam,options)
    return nothing
end

Base.@specialize
include("database_utils.jl")

"""
    params, sites = getparams(components,locations;kwargs...)

returns a `Dict{String,ClapeyronParam}` containing all the parameters found for the list of components
in the available CSVs. `locations` are the locations relative to `Clapeyron` database. the available keywords are the ones used in [`ParamOptions`](@ref)

if `return_sites` is set to false, `getparams` will only return the found params.
"""
function getparams(components, 
                    locations::Array{String,1}=String[]; 
                    userlocations::Vector{String}=String[],
                    asymmetricparams::Vector{String}=String[],
                    ignore_missing_singleparams::Vector{String}=String[],
                    ignore_headers::Vector{String} =  ["dipprnumber", "smiles"],
                    verbose::Bool=false,
                    species_columnreference::String="species",
                    source_columnreference::String="source",
                    site_columnreference::String="site",
                    group_columnreference::String="groups",
                    normalisecomponents::Bool=true,
                    return_sites::Bool = true,
                    component_delimiter::String = "~|~"
                    )

    options = ParamOptions(;userlocations,
                            asymmetricparams,
                            ignore_missing_singleparams,
                            ignore_headers,
                            verbose,
                            species_columnreference,
                            source_columnreference,
                            site_columnreference,
                            group_columnreference,
                            normalisecomponents,
                            return_sites,
                            component_delimiter)
    
    # locations is a list of paths relative to the Clapeyron database directory.
    # userlocations is a list of paths input by the user.
    # If parameters exist in multiple files, Clapeyron gives priority to files in later paths.
    # asymmetricparams is a list of parameters for which matrix reflection is disabled.
    # ignore_missingsingleparams gives users the option to disable component existence check in single params.                   
    return getparams(components,locations,options)
end
function getparams(components::Vector{String},locations::Vector{String},options::ParamOptions)
    #generate one string of params
    filepaths = flattenfilepaths(locations,options.userlocations)
    #merge all found params
    allparams,allnotfoundparams = createparams(components, filepaths, options)
    #generate ClapeyronParams
    result, allcomponentsites = compile_params(components,allparams,allnotfoundparams,options)
    #check values
    for v in values(result) 
        is_valid_param(v,options)
    end

    if !(options.return_sites)
        return result
    end
    if any(x isa AssocParam for x in values(result))
        sites = buildsites(result,components,allcomponentsites,options)
        return result,sites
    else
        return result
    end
end

function buildsites(result,components,allcomponentsites,options)
    n_sites_columns = options.n_sites_columns
    v = String[]
    for sitei in allcomponentsites
        append!(v,sitei)
    end
    unique!(v)
    iszero(length(v)) && return SiteParam(components)
    if !any(haskey(result,n_sites_columns[vi]) for vi in v)
        @error """No columns containing number of sites were found. Supposing zero sites.
        If your model doesn't use sites, but some input CSV folders contain Association params. consider using return_sites=false.
        If there are columns containing number of sites, then those aren't recognized.
        Consider passing a Dict with the mappings between the sites and the name of the column containing the number of said sites,
        with the n_sites_columns keyword"
        """
        assoc_csv = Set(String[])
        for x in values(result)
            if x isa AssocParam
                for csvx in x.sourcecsvs
                push!(assoc_csv,csvx)
                end
            end
        end
        assoc_csv = collect(assoc_csv)
        @error "Parsed Association CSV were:"
        for csv in assoc_csv
            println(csv)
        end
        return SiteParam(components)
    end
    
    n_sites_dict = Dict{String,SingleParam{Int}}(vi => result[n_sites_columns[vi]] for vi in v)    
    return SiteParam(n_sites_dict,allcomponentsites)
end

function getparams(groups::GroupParam, locations::Vector{String}=String[],options::ParamOptions=DefaultOptions)
    return getparams(groups.flattenedgroups, locations,options)
end

function getparams(components::String, locations::Vector{String}=String[],options::ParamOptions=DefaultOptions)
    return getparams([components],locations,options)
end

function findsites(data::Dict,components::Vector;verbose = false)
    sites = Dict(components .=> [Set{String}() for _ ∈ 1:length(components)])
    for raw in values(data)
        if raw.type === assocdata
            for (c1,c2,s1,s2) in raw.component_info
            push!(sites[c1], s1)
            push!(sites[c2], s2)
            end
        end
    end
    output = Array{Array{String,1}}(undef, 0)
    for component ∈ components
        push!(output, collect(sites[component]))
    end
    verbose && @info("Found sites for $components are $(output).")
    return output
end

function createparams(components::Vector{String},
                    filepaths::Vector{String},
                    options::ParamOptions = DefaultOptions,
                    parsegroups = false)

    allparams = Dict{String,RawParam}()
    allnotfoundparams = Dict{String,CSVType}()
    # Read the filepaths in reverse in order to ensure that unused sources do not get added.
    for filepath ∈ reverse(filepaths)
        csvtype = readcsvtype(filepath)
        if csvtype == groupdata && !parsegroups
            continue
        end
        foundparams, notfoundparams = findparamsincsv(components, filepath,options,parsegroups)
        
        #Merge found data
        for vv in foundparams
            kk = vv.name
            if haskey(allparams,kk)
                vv2 = allparams[kk]
                vv = joindata!(vv2,vv)
            end
            allparams[kk] = vv
        end
        #Merge not found data
        for (kk,vv) in pairs(notfoundparams)
            if haskey(allnotfoundparams,kk)
                vv2 = allnotfoundparams[kk]
                vv, success = joindata!(vv2,vv)
                !success && error_clashing_headers(vv2,vv,kk) #Clashing headers error
            end
            allnotfoundparams[kk] = vv
        end
    end
    #delete all found params from allnotfoundparams
    for (kk,vv) in allparams
            delete!(allnotfoundparams,kk)
    end

    return allparams,allnotfoundparams
end

function compile_params(components,allparams,allnotfoundparams,options)
    #Generate component sites with the RawParam tapes
    allcomponentsites = findsites(allparams,components)
    #Compile Params
    result = Dict{String,ClapeyronParam}(k => compile_param(components,k,v,allcomponentsites,options) for (k,v) in allparams)
    for (kk,vv) in allnotfoundparams
        result[kk] = compile_param(components,kk,vv,allcomponentsites,options) 
    end
    return result,allcomponentsites
end

@noinline function _col_indices_error(header)
    throw(error("Header ", normalised_columnreference, " not found."))
end

function col_indices(csvtype,headernames,options=DefaultOptions)
    columnreference = options.species_columnreference
    normalised_columnreference = normalisestring(columnreference)

    idx_species = 0
    idx_groups = 0
    idx_species1 = 0
    idx_species2 = 0
    idx_sites1 = 0
    idx_sites2 = 0

    if csvtype === singledata || csvtype === groupdata
        lookupcolumnindex = findfirst(isequal(normalised_columnreference), headernames)
        isnothing(lookupcolumnindex) && _col_indices_error(normalised_columnreference)
        idx_species = lookupcolumnindex
        if csvtype === groupdata
            groupcolumnreference = options.group_columnreference
            normalised_groupcolumnreference = normalisestring(groupcolumnreference)
            lookupgroupcolumnindex = findfirst(isequal(normalised_groupcolumnreference), headernames)
            isnothing(lookupgroupcolumnindex) && _col_indices_error(normalised_groupcolumnreference)
            idx_groups = lookupgroupcolumnindex
        end
    
    elseif csvtype === pairdata || csvtype == assocdata
        normalised_columnreference1 = normalised_columnreference * '1'
        normalised_columnreference2 = normalised_columnreference * '2'
        lookupcolumnindex1 = findfirst(isequal(normalised_columnreference1), headernames)
        lookupcolumnindex2 = findfirst(isequal(normalised_columnreference2), headernames)
        isnothing(lookupcolumnindex1) && _col_indices_error(normalised_columnreference1)
        isnothing(lookupcolumnindex2) && _col_indices_error(normalised_columnreference2)
        idx_species1 = lookupcolumnindex1
        idx_species2 = lookupcolumnindex2
        if csvtype == assocdata
            sitecolumnreference = options.site_columnreference
            normalised_sitecolumnreference = normalisestring(sitecolumnreference)
            normalised_sitecolumnreference1 = normalised_sitecolumnreference * '1'
            normalised_sitecolumnreference2 = normalised_sitecolumnreference * '2'
            lookupsitecolumnindex1 = findfirst(isequal(normalised_sitecolumnreference1), headernames)
            lookupsitecolumnindex2 = findfirst(isequal(normalised_sitecolumnreference2), headernames)
            isnothing(lookupsitecolumnindex1) && _col_indices_error(normalised_sitecolumnreference1)
            isnothing(lookupsitecolumnindex2) && _col_indices_error(normalised_sitecolumnreference2)
            idx_sites1 = lookupsitecolumnindex1
            idx_sites2 = lookupsitecolumnindex2
        end
    end

    _single = (idx_species,idx_groups)
    _pair = (idx_species1,idx_species2)
    _assoc = (idx_sites1,idx_sites2)
    return (_single,_pair,_assoc)
end

function findparamsincsv(components,filepath,options::ParamOptions = DefaultOptions,parsegroups = false)

    headerparams = readheaderparams(filepath,options)
    sourcecolumnreference = options.source_columnreference
    verbose = options.verbose
    normalisecomponents = options.normalisecomponents
    component_delimiter = options.component_delimiter

    csvtype = readcsvtype(filepath)
    grouptype = "" #TODO: actually parse this from csvs
    df = CSV.File(filepath; header=3, pool=0,silencewarnings=true)
    csvheaders = String.(Tables.columnnames(df))
    normalised_components = normalisestring.(components,normalisecomponents)
    components_dict = Dict(v => k for (k,v) in pairs(normalised_components))
    normalised_csvheaders = normalisestring.(csvheaders)
    normalised_headerparams = normalisestring.(headerparams)
    if normalised_headerparams ⊈ normalised_csvheaders 
        error("Headers ", setdiff(normalised_headerparams, normalised_csvheaders), " not present in csv header.")
    end

    #function output
    foundvalues = Vector{RawParam}(undef,0)
    notfoundvalues = Dict{String,CSVType}(headerparam => csvtype for headerparam in headerparams)

    normalised_sourcecolumnreference = normalisestring(sourcecolumnreference)
    getsources = false
    if normalised_sourcecolumnreference ∈ normalised_csvheaders
        getsources = true
        sourcecolumn = findfirst(isequal(normalised_sourcecolumnreference), normalised_csvheaders)
    end

    single_idx,pair_idx,assoc_idx = col_indices(csvtype,normalised_csvheaders,options)
    lookupcolumnindex,groupindex = single_idx
    lookupcolumnindex1,lookupcolumnindex2 = pair_idx
    lookupsitecolumnindex1,lookupsitecolumnindex2 = assoc_idx
    headerparams_indices = [findfirst(isequal(i),normalised_csvheaders) for i in normalised_headerparams]
    lookupcolumnindex = max(lookupcolumnindex,lookupcolumnindex1)
    
    verbose && __verbose_findparams_start(filepath,components,headerparams,parsegroups,csvtype)
    #list of all species
    species_list = normalisestring.(Tables.getcolumn(df,lookupcolumnindex),normalisecomponents)
    
    #indices where data could be (they could be missing)
    #on pair and assoc, this is just the first component, we need to reduce the valid indices again
    found_indices0,comp_indices = _indexin(components_dict,species_list,component_delimiter,1:length(species_list))
    dfR = Tables.rows(df)
    EMPTY_STR = ""

    if csvtype == singledata || ((csvtype == groupdata) && parsegroups)
        found_indices = found_indices0
        l = length(found_indices)
        if l != 0
            _data = dfR[found_indices]
            _comp = [(components[c],EMPTY_STR,EMPTY_STR,EMPTY_STR) for c in comp_indices]
            _sources = fill(EMPTY_STR,l)
            _csv = fill(String(filepath),l)
        end
    elseif csvtype == pairdata && !parsegroups
        species2_list = normalisestring.(Tables.getcolumn(df,lookupcolumnindex2)[found_indices0],normalisecomponents)
        found_indices2,comp_indices2 = _indexin(components_dict,species2_list,component_delimiter,1:length(species2_list))
        comp_indices1 = comp_indices[found_indices2]
        found_indices2 = found_indices0[found_indices2]
        l = length(found_indices2)
        if l != 0
            _data = dfR[found_indices2]
            _comp = [(components[c1],components[c2],EMPTY_STR,EMPTY_STR) for (c1,c2) in zip(comp_indices1,comp_indices2)] 
            _sources = fill(EMPTY_STR,l)
            _csv = fill(String(filepath),l)
        end
    elseif csvtype == assocdata && !parsegroups  
        species2_list = normalisestring.(Tables.getcolumn(df,lookupcolumnindex2)[found_indices0],normalisecomponents)        
        found_indices2,comp_indices2 = _indexin(components_dict,species2_list,component_delimiter,1:length(species2_list))
        comp_indices1 = comp_indices[found_indices2]
        found_indices2 = found_indices0[found_indices2]
        l = length(found_indices2)
        if l != 0
            _data = dfR[found_indices2]
            _site1 = getindex.(_data,lookupsitecolumnindex1)
            _site2 = getindex.(_data,lookupsitecolumnindex2)
            _comp = [(components[c1],components[c2],String(s1),String(s2)) for (c1,c2,s1,s2) in zip(comp_indices1,comp_indices2,_site1,_site2)] 
            _sources = fill(EMPTY_STR,l)
            _csv = fill(String(filepath),l)
        end
    elseif csvtype == groupdata && !parsegroups 
        return foundvalues, notfoundvalues
    else
        error("File is of type ", String(csvtype), " and cannot be read with this function.")
    end

    if l != 0
        #if getsources, then we actually put the sources in inside the preallocated _sources vector
        if getsources
            _raw_sources = getindex.(_data,sourcecolumn)
            for i in eachindex(_sources)
                source_i = _raw_sources[i]
                !ismissing(source_i) && (_sources[i] = source_i) 
            end
        end
        #with the raw data preallocated, we now store it in a RawParam.
        for (headerparam,idx) ∈ zip(headerparams,headerparams_indices)
            _vals = getindex.(_data,idx)
            s = findall(!ismissing,_vals) #filter nonmissing values
            if !iszero(s)
                __vals = [_vals[i] for i in s] #removes the missing type and eliminates bitvectors
                __sources = [_sources[i] for i in s]
                __csv = [_csv[i] for i in s]
                raw = RawParam(headerparam,_comp[s],__vals,__sources,__csv,csvtype,grouptype)
                push!(foundvalues,raw)
            end
        end
    end
    #store all headers that didn't had a result.
    for rawparam in foundvalues
        delete!(notfoundvalues,rawparam.name)
    end
    
    verbose && __verbose_findparams_found(foundvalues) #print all found values

    return foundvalues, notfoundvalues
end

#verbose functionality, is executed for each csv when verbose == true

function __assoc_string(pair)
    "($(pair[1]),$(pair[3])) ⇋ ($(pair[2]), $(pair[4]))"
end

function __verbose_findparams_start(filepath,components,headerparams,parsegroups,csvtype)
    csv_string = Symbol(csvtype)
    if !parsegroups
        if csvtype == groupdata
            @info("Skipping $csv_string csv $filepath")    
        else
            @info("Searching for $csv_string headers $headerparams for query $components at $filepath ...")
        end
    else
        if csvtype == groupdata 
            @info("Searching for groups for components $components at $filepath ...")
        else
            @info("Skipping $csv_string csv $filepath")  
        end
    end
end


function __verbose_findparams_found(foundvalues)
    c = Base.text_colors
    reset = c[:normal]
    blue = c[:bold] * c[:cyan]
    for (k,v) in pairs(foundvalues)
        if v.type == singledata
            vdict = Dict(pair[1] => val for (pair,val) in zip(v.component_info,v.data))
            kk = blue * v.name * reset
            @info("""Found single component data: $kk with values: 
            $vdict
            """)
        elseif v.type == pairdata
            vdict = Dict((pair[1],pair[2]) => val for (pair,val) in zip(v.component_info,v.data))
            kk = blue * v.name * reset
            @info("""Found pair component data: $kk with values: 
            $vdict
            """)
        elseif v.type == assocdata
            vdict = Dict(__assoc_string(pair) => val for (pair,val) in zip(v.component_info,v.data))
            kk = blue * v.name * reset
            @info("""Found association component data: $kk with values: 
            $vdict
            """)
        elseif v.type == groupdata
            vdict = Dict(pair[1] => val for (pair,val) in zip(v.component_info,v.data))
            #println(val) for val in v.data)
            @info("""Found group data: 
            $vdict
            """)
        end
    end
end

const readcsvtype_keywords  = ["like", "single", "unlike", "pair", "assoc", "group", "groups"]

function readcsvtype(filepath)
    # Searches for type from second line of CSV.
    keywords = readcsvtype_keywords
    words = split(lowercase(rstrip(getline(String(filepath), 2), ',')), ' ')
    foundkeywords = intersect(words, keywords)
    if isempty(foundkeywords) 
        error("Unable to determine type of database", filepath, ". Check that keyword is present on Line 2.")
    end
    if length(foundkeywords) > 1 
        error("Multiple keywords found in database ", filepath, ": ", foundkeywords)
    end
    _readcsvtype(only(foundkeywords)) 
end

function _readcsvtype(key)
    key == "single" && return singledata
    key == "like" && return singledata
    key == "pair" && return pairdata
    key == "unlike" && return pairdata
    key == "assoc" && return assocdata
    key == "group" && return groupdata
    key == "groups" && return groupdata
    error("Unable to determine database type of $key")
end

function readheaderparams(filepath::AbstractString, options::ParamOptions = DefaultOptions,headerline::Int = 3)
    # Returns array of filtered header strings at line 3.
    ignorelist = deepcopy(options.ignore_headers)
    push!(ignorelist,options.species_columnreference)
    push!(ignorelist,options.source_columnreference)
    push!(ignorelist,options.site_columnreference)
    headers = split(getline(filepath, headerline), ',',keepempty = false)
    while last(headers) == ""
        pop!(headers)
    end
    return String.(filter(x -> normalisestring(x; tofilter=r"[ \-\_\d]") ∉ ignorelist, headers))
end

"""
    singletopair(params::Vector,outputmissing=zero(T))

Generates a square matrix, filled with "zeros" (considering the "zero" of a string, a empty string). 
The generated matrix will have the values of `params` in the diagonal.
If missing is passed, the matrix will be filled with `missing`
"""
function singletopair(params::Vector{T1},::T2 =_zero(T1)) where {T1,T2}
    len = length(params)
    T = Union{T1,T2}
    output = Matrix{T}(undef,len,len)
    fill!(output,_zero(T))
    @inbounds  for i in 1:len
        output[i,i] = params[i]
    end
    return output
end

function GroupParam(gccomponents::Vector, 
    grouplocations::Vector{String}=String[]; 
    usergrouplocations::Vector{String}=String[], 
    verbose::Bool=false)
    options = ParamOptions(;usergrouplocations,verbose)
    return GroupParam(gccomponents,grouplocations,options)        
end 

function GroupParam(gccomponents, 
    grouplocations::Array{String,1}=String[],
    options::ParamOptions = DefaultOptions)
    # The format for gccomponents is an arary of either the species name (if it
    # available in the Clapeyron database, or a tuple consisting of the species
    # name, followed by a list of group => multiplicity pairs.  For example:
    # gccomponents = ["ethane",
    #                ("hexane", ["CH3" => 2, "CH2" => 4]),
    #                ("octane", ["CH3" => 2, "CH2" => 6])]
    components = Vector{String}(undef,length(gccomponents))
    
    to_lookup = fill(false,length(components))
    for (i,gcpair) in pairs(gccomponents)
        if gcpair isa Tuple{String}
            components[i]  = only(gcpair)
            to_lookup[i] = true
        elseif gcpair isa String
            components[i]  = gcpair
            to_lookup[i] = true
        elseif gcpair isa Tuple{String, Vector{Pair{String, Int64}}}
            components[i]  = first(pair)
        else
            error("The format of the components is incorrect.")
        end
    end
    #using parsing machinery
    usergrouplocations = options.usergrouplocations
    componentstolookup = components[to_lookup]
    filepaths = flattenfilepaths(grouplocations,usergrouplocations) 
    allparams,allnotfoundparams = createparams(componentstolookup, filepaths, options, true) #merge all found params
    raw_result, _ = compile_params(componentstolookup,allparams,allnotfoundparams,options) #generate ClapeyronParams
    raw_groups = raw_result["groups"] #SingleParam{String}
    is_valid_param(raw_groups,options) #this will check if we actually found all params, via single missing detection.
    groupsourcecsvs = raw_groups.sourcecsvs
  
    gccomponents_parsed = PARSED_GROUP_VECTOR_TYPE(undef,length(gccomponents))
    j = 0
    for (i,needs_to_parse_group_i) ∈ pairs(to_lookup)
        if needs_to_parse_group_i #we looked up this component, and if we are here, it exists.
            j += 1
            gcdata = _parse_group_string(raw_groups.values[j])
            gccomponents_parsed[i] = (components[i],gcdata)
        else
            gccomponents_parsed[i] = gccomponents[i]
        end
    end
    return GroupParam(gccomponents_parsed,groupsourcecsvs)
end

function _parse_group_string(gc::String)
    gc_strip = strip(gc)
    if startswith(gc_strip,"[") && endswith(gc_strip,"]")
        gc_without_brackets = chop(gc_strip,head = 1,tail = 1)
        gcpairs = split(gc_without_brackets,",")
        result = Vector{Pair{String,Int}}(undef,length(gcpairs))
        for (i,gcpair) in pairs(gcpairs)
            x = split(strip(gcpair),"=>") # """ "__group__"__=>__number__""" 
            length(x) != 2 && throw(error("incorrect group format"))
            raw_group_i = strip(x[1])
            if startswith(raw_group_i,"\"") && endswith(raw_group_i,"\"")
                group_i = chop(raw_group_i,head = 1,tail = 1)
                num = parse(Int64,x[2])
                result[i] = group_i => num
            else
                throw(error("incorrect group format"))
            end
        end
        return result
    else
        throw(error("incorrect group format"))
    end

end
export getparams