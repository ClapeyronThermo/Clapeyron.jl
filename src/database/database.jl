@enum CSVType invaliddata singledata pairdata assocdata groupdata

include("database_rawparam.jl")
include("database_utils.jl")

"""
    params, sites = getparams(components,locations;kwargs...)

returns a `Dict{String,ClapeyronParam}` containing all the parameters found for the list of components
in the available CSVs. `locations` are the locations relative to `Clapeyron` database. the available keywords are the ones used ∈ [`ParamOptions`](@ref)

if `return_sites` is set to false, `getparams` will only return the found params.
"""
function getparams(components,
                    locations::Array{String,1}=String[];
                    userlocations::Vector{String}=String[],
                    asymmetricparams::Vector{String}=String[],
                    ignore_missing_singleparams::Vector{String}=String[],
                    ignore_headers::Vector{String} =  ["dipprnumber", "smiles", "cas"],
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
    # If parameters exist ∈ multiple files, Clapeyron gives priority to files ∈ later paths.
    # asymmetricparams is a list of parameters for which matrix reflection is disabled.
    # ignore_missingsingleparams gives users the option to disable component existence check ∈ single params.
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
    for v ∈ values(result)
        is_valid_param(v,options)
    end

    if !(options.return_sites)
        return result
    end
    if any(x isa AssocParam for x ∈ values(result))
        sites = buildsites(result,components,allcomponentsites,options)
        return result,sites
    else
        return result
    end
end

function buildsites(result,components,allcomponentsites,options)
    n_sites_columns = options.n_sites_columns
    v = String[]
    for sitei ∈ allcomponentsites
        append!(v,sitei)
    end
    unique!(v)
    iszero(length(v)) && return SiteParam(components)
    if !any(haskey(result,n_sites_columns[vi]) for vi ∈ v)
       return __warning_no_site_vals(result,components)
    end

    n_sites_dict = Dict{String,SingleParam{Int}}(vi => result[n_sites_columns[vi]] for vi ∈ v)
    return SiteParam(n_sites_dict,allcomponentsites)
end

function __warning_no_site_vals(result,components)
    @error """No columns containing number of sites were found. Supposing zero sites.
    If your model doesn't use sites, but some input CSV folders contain Association params. consider using return_sites=false.
    If there are columns containing number of sites, then those aren't recognized.
    Consider passing a Dict with the mappings between the sites and the name of the column containing the number of said sites,
    with the n_sites_columns keyword"
    """
    assoc_csv = Set(String[])
    for x ∈ values(result)
        if x isa AssocParam
            for csvx ∈ x.sourcecsvs
            push!(assoc_csv,csvx)
            end
        end
    end
    assoc_csv = collect(assoc_csv)
    @error "Parsed Association CSV were:"
    for csv ∈ assoc_csv
        println(csv)
    end
    return SiteParam(components)
end

function getparams(groups::GroupParam, locations::Vector{String}=String[],options::ParamOptions=DefaultOptions)
    return getparams(groups.flattenedgroups, locations,options)
end

function getparams(components::String, locations::Vector{String}=String[],options::ParamOptions=DefaultOptions)
    return getparams([components],locations,options)
end

function findsites(data::Dict,components::Vector;verbose = false)
    sites = Dict(components .=> [Set{String}() for _ ∈ 1:length(components)])
    for raw ∈ values(data)
        if raw.type === assocdata
            for (c1,c2,s1,s2) ∈ raw.component_info
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
@nospecialize
function createparams(components::Vector{String},
                    filepaths::Vector{String},
                    options::ParamOptions = DefaultOptions,
                    parsegroups = false)

    allparams = Dict{String,RawParam}()
    allnotfoundparams = Dict{String,CSVType}()
    for filepath ∈ filepaths
        
        _replace = startswith(filepath,"@REPLACE")
        if _replace
            filepath = chop(filepath,head = 9, tail = 0)
        end
        
        csv_options = read_csv_options(filepath)
        csvtype = csv_options.csvtype
        if csvtype == groupdata && !parsegroups
            continue
        end
        if csvtype == invaliddata
            if options.verbose
                __verbose_findparams_invaliddata(filepath)
            end
            continue
        end
        foundparams, notfoundparams = findparamsincsv(components,filepath,options,parsegroups,csv_options)

        #Merge found data
        for vv ∈ foundparams
            kk = vv.name
            #we merge if the filepath is not set to replace the current values
            if haskey(allparams,kk) && !_replace
                vv2 = allparams[kk]
                vv = joindata!(vv2,vv)
            end
            allparams[kk] = vv
        end
        #Merge not found data
        for (kk,vv) ∈ pairs(notfoundparams)
            if haskey(allnotfoundparams,kk)
                vv2 = allnotfoundparams[kk]
                vv, success = joindata!(vv2,vv)
                !success && error_clashing_headers(vv2,vv,kk) #Clashing headers error
            end
            allnotfoundparams[kk] = vv
        end

        if _replace #if the paramter is not found, that means that we want to erase that param.
            for (kk,vv) ∈ pairs(notfoundparams)
                delete!(allparams,kk)
            end
        end
    end
    #delete all found params from allnotfoundparams
    for (kk,vv) ∈ allparams
        delete!(allnotfoundparams,kk)
    end

    return allparams,allnotfoundparams
end
@specialize
function compile_params(components,allparams,allnotfoundparams,options)
    #Generate component sites with the RawParam tapes
    allcomponentsites = findsites(allparams,components)
    #Compile Params
    result = Dict{String,ClapeyronParam}(k => compile_param(components,k,v,allcomponentsites,options) for (k,v) ∈ allparams)
    for (kk,vv) ∈ allnotfoundparams
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

function findparamsincsv(components,filepath,options::ParamOptions = DefaultOptions,parsegroups = false,csv_file_options = read_csv_options(filepath))

    headerparams = readheaderparams(filepath,options)
    sourcecolumnreference = options.source_columnreference
    verbose = options.verbose
    normalisecomponents = options.normalisecomponents
    component_delimiter = options.component_delimiter

    grouptype = csv_file_options.grouptype
    csvtype = csv_file_options.csvtype
    if startswith(filepath,"Clapeyron Database File")
        df = CSV.File(IOBuffer(filepath); header=3, pool=0,silencewarnings=true)
    else
        df = CSV.File(filepath; header=3, pool=0,silencewarnings=true)

    end
    csvheaders = String.(Tables.columnnames(df))
    normalised_components = normalisestring.(components,normalisecomponents)
    components_dict = Dict(v => k for (k,v) ∈ pairs(normalised_components))
    normalised_csvheaders = normalisestring.(csvheaders)
    normalised_headerparams = normalisestring.(headerparams)
    if normalised_headerparams ⊈ normalised_csvheaders
        error("Headers ", setdiff(normalised_headerparams, normalised_csvheaders), " not present ∈ csv header.")
    end

    #function output
    foundvalues = Vector{RawParam}(undef,0)
    notfoundvalues = Dict{String,CSVType}(headerparam => csvtype for headerparam ∈ headerparams)

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
    headerparams_indices = [findfirst(isequal(i),normalised_csvheaders) for i ∈ normalised_headerparams]
    lookupcolumnindex = max(lookupcolumnindex,lookupcolumnindex1)

    verbose && __verbose_findparams_start(filepath,components,headerparams,parsegroups,csvtype)
    #list of all species
    species_list = normalisestring.(Tables.getcolumn(df,lookupcolumnindex),normalisecomponents)

    #indices where data could be (they could be missing)
    #on pair and assoc, this is just the first component, we need to reduce the valid indices again
    found_indices0,comp_indices = _indexin(components_dict,species_list,component_delimiter,1:length(species_list))
    dfR = Tables.rowtable(df)
    EMPTY_STR = ""

    if csvtype == singledata || ((csvtype == groupdata) && parsegroups)
        found_indices = found_indices0
        l = length(found_indices)
        if l != 0
            _data = dfR[found_indices]
            _comp = [(components[c],EMPTY_STR,EMPTY_STR,EMPTY_STR) for c ∈ comp_indices]
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
            _comp = [(components[c1],components[c2],EMPTY_STR,EMPTY_STR) for (c1,c2) ∈ zip(comp_indices1,comp_indices2)]
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
            _site1 = identity.(getindex.(_data,lookupsitecolumnindex1))
            _site2 = identity.(getindex.(_data,lookupsitecolumnindex2))
            _comp = [(components[c1],components[c2],String(s1),String(s2)) for (c1,c2,s1,s2) ∈ zip(comp_indices1,comp_indices2,_site1,_site2)]
            _sources = fill(EMPTY_STR,l)
            _csv = fill(String(filepath),l)
        end
    elseif csvtype == groupdata && !parsegroups
        return foundvalues, notfoundvalues
    else
        error("File is of type ", String(csvtype), " and cannot be read with this function.")
    end

    if l != 0
        #if getsources, then we actually put the sources ∈ inside the preallocated _sources vector
        if getsources
            _raw_sources = getindex.(_data,sourcecolumn)
            for i ∈ eachindex(_sources)
                source_i = _raw_sources[i]
                !ismissing(source_i) && (_sources[i] = source_i)
            end
        end
        #with the raw data preallocated, we now store it ∈ a RawParam.
        for (headerparam,idx) ∈ zip(headerparams,headerparams_indices)
            _vals = getindex.(_data,idx)
            s = findall(!ismissing,_vals) #filter nonmissing values
            if !iszero(length(s))
                __vals = [_vals[i] for i ∈ s] #removes the missing type and eliminates bitvectors
                __sources = [_sources[i] for i ∈ s]
                __csv = [_csv[i] for i ∈ s]
                raw = RawParam(headerparam,_comp[s],__vals,__sources,__csv,csvtype,grouptype)
                push!(foundvalues,raw)
            end
        end
    end
    #store all headers that didn't had a result.
    for rawparam ∈ foundvalues
        delete!(notfoundvalues,rawparam.name)
    end

    verbose && __verbose_findparams_found(foundvalues) #print all found values

    return foundvalues, notfoundvalues
end

#verbose functionality, is executed for each csv when verbose == true

function __verbose_findparams_invaliddata(filepath)
    @warn "Skipping $filepath, cannot infer correct csv type. Check line 2 of the CSV to see if it has valid information."
end

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
    for v ∈ foundvalues
        if v.type == singledata
            vdict = Dict(pair[1] => val for (pair,val) ∈ zip(v.component_info,v.data))
            kk = info_color(v.name)
            @info("""Found single component data: $kk with values:
            $vdict
            """)
        elseif v.type == pairdata
            vdict = Dict((pair[1],pair[2]) => val for (pair,val) ∈ zip(v.component_info,v.data))
            kk = info_color(v.name)
            @info("""Found pair component data: $kk with values:
            $vdict
            """)
        elseif v.type == assocdata
            vdict = Dict(__assoc_string(pair) => val for (pair,val) ∈ zip(v.component_info,v.data))
            kk = info_color(v.name)
            @info("""Found association component data: $kk with values:
            $vdict
            """)
        elseif v.type == groupdata
            #println(val) for val ∈ v.data)
            vals = "Dict("
            for (pair,val) ∈ zip(v.component_info,v.data)
                pairi = pair[1] * " => " * string(val)
                vals = vals * pairi * ", "
            end
            vals = chop(vals,tail=2) *")"
            @info("""Found group data:
            $vals
            """)
        end
    end
end

const readcsvtype_keywords  = ["like", "single", "unlike", "pair", "assoc", "association", "group", "groups"]

function read_csv_options(filepath)
    line = getline(String(filepath), 2)
    re = r"\[.*\]"
    maybe_opts = match(re,line)

    # Searches for type from second line of CSV.
    has_csv_options = !isnothing(maybe_opts)
    if has_csv_options
        opts = chop(maybe_opts.match,head = 1,tail = 1)
        return __get_options(opts)
    else
        keywords = readcsvtype_keywords
        words = split(lowercase(strip(line, ',')), ' ')
        foundkeywords = intersect(words, keywords)
        return (csvtype = _readcsvtype(foundkeywords),grouptype = :unknown)
    end
end


function _readcsvtype(collection)
    length(collection) != 1 && return invaliddata
    key = only(collection)
    return _readcsvtype(key)
end
function _readcsvtype(key::AbstractString)
    key == "single" && return singledata
    key == "like" && return singledata
    key == "pair" && return pairdata
    key == "unlike" && return pairdata
    key == "assoc" && return assocdata
    key == "group" && return groupdata
    key == "groups" && return groupdata
    key == "invalid" && return invaliddata
    return invaliddata
end

function __get_options(data)
    opts = split(data,',')
    opts_dict = Dict{String,String}()
    for opt in opts
        k,v = split(replace(opt,' ' => ""),'=')
        opts_dict[k] = v
    end
    _csvtype = _readcsvtype(get(opts_dict,"csvtype","invalid"))
    _grouptype = Symbol(get(opts_dict,"grouptype","unkwown"))
    return (csvtype = _csvtype,grouptype = _grouptype)
end

function readheaderparams(filepath::AbstractString, options::ParamOptions = DefaultOptions,headerline::Int = 3)
    # Returns array of filtered header strings at line 3.
    ignorelist = deepcopy(options.ignore_headers)
    push!(ignorelist,options.species_columnreference)
    push!(ignorelist,options.source_columnreference)
    push!(ignorelist,options.site_columnreference)
    headertext = getline(filepath, headerline)
    headers = split(rstrip(headertext,','), ',',keepempty = false)
    return String.(filter(x -> normalisestring(x; tofilter=r"[ \-\_\d]") ∉ ignorelist, headers))
end

function GroupParam(gccomponents::Vector,
    group_locations::Vector{String}=String[];
    group_userlocations::Vector{String}=String[],
    verbose::Bool=false)
    options = ParamOptions(;group_userlocations,verbose)
    return GroupParam(gccomponents,group_locations,options)
end

function GroupParam(gccomponents,
    grouplocations::Array{String,1}=String[],
    options::ParamOptions = DefaultOptions)
    # The format for gccomponents is an arary of either the species name (if it
    # available ∈ the Clapeyron database, or a tuple consisting of the species
    # name, followed by a list of group => multiplicity pairs.  For example:
    # gccomponents = ["ethane",
    #                ("hexane", ["CH3" => 2, "CH2" => 4]),
    #                ("octane", ["CH3" => 2, "CH2" => 6])]
    components = Vector{String}(undef,length(gccomponents))

    to_lookup = fill(false,length(components))
    for (i,gcpair) ∈ pairs(gccomponents)
        if gcpair isa Tuple{String}
            components[i]  = only(gcpair)
            to_lookup[i] = true
        elseif gcpair isa String
            components[i]  = gcpair
            to_lookup[i] = true
        elseif gcpair isa Tuple{String, Vector{Pair{String, Int64}}}
            components[i]  = first(gcpair)
        else
            error("The format of the components is incorrect.")
        end
    end
    #using parsing machinery
    usergrouplocations = options.group_userlocations
    componentstolookup = components[to_lookup]
    filepaths = flattenfilepaths(grouplocations,usergrouplocations)
    allparams,allnotfoundparams = createparams(componentstolookup, filepaths, options, true) #merge all found params
    raw_result, _ = compile_params(componentstolookup,allparams,allnotfoundparams,options) #generate ClapeyronParams
    raw_groups = raw_result["groups"] #SingleParam{String}
    is_valid_param(raw_groups,options) #this will check if we actually found all params, via single missing detection.
    groupsourcecsvs = raw_groups.sourcecsvs
    
    ##TODO: extract groups
    if haskey(allparams,"groups")
        grouptype = allparams["groups"].grouptype
    else
        grouptype = :unknown
    end

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
    return GroupParam(gccomponents_parsed,grouptype,groupsourcecsvs)
end

function _parse_group_string(gc::String)
    gc_strip = strip(gc)
    if startswith(gc_strip,"[") && endswith(gc_strip,"]")
        gc_without_brackets = chop(gc_strip,head = 1,tail = 1)
        gcpairs = split(gc_without_brackets,",")
        result = Vector{Pair{String,Int}}(undef,length(gcpairs))
        for (i,gcpair) ∈ pairs(gcpairs)
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