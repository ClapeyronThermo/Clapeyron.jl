function fast_parse_grouptype(filepaths::Vector{String})
    #only parses grouptype, if present in any CSV, is used. if not, return unkwown
    grouptype = :not_set
   
    for filepath ∈ filepaths
        _replace = startswith(filepath,"@REPLACE")
        if _replace
            filepath = chop(filepath,head = 9, tail = 0)
        end
        csv_options = read_csv_options(filepath)
        new_grouptype = csv_options.grouptype
        if grouptype == :not_set
            grouptype = new_grouptype
        else
            if grouptype !== new_grouptype
                if new_grouptype != :unkwown #for backwards compatibility
                    error_different_grouptype(grouptype,new_grouptype)
                end
            end
        end
    end
    grouptype == :not_set && (grouptype = :unkwown)
    return grouptype

end

#used to parse """["CH => 2, "OH" = 2, "CH3" => 2]"""
function _parse_group_string(gc::String,gctype=String)
    gc_strip = strip(gc)
    if startswith(gc_strip,"[") && endswith(gc_strip,"]")
        gc_without_brackets = chop(gc_strip,head = 1,tail = 1)
        gcpairs = split(gc_without_brackets,",")
        result = Vector{Pair{gctype,Int}}(undef,length(gcpairs))
        for (i,gcpair) ∈ pairs(gcpairs)
            raw_group_i,raw_num = _parse_kv(gcpair,"=>")
            if startswith(raw_group_i,"\"") && endswith(raw_group_i,"\"")
                group_i = _parse_group_string_key(raw_group_i,gctype)
                num = parse(Int64,raw_num)
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

#"CH3"
function _parse_group_string_key(k,::Type{String})
    return chop(k,head = 1,tail = 1)
end

#("CH3","CH2")
function _parse_group_string_key(k,::Type{NTuple{2,String}})
     kk = chop(strip(k),head = 1,tail = 1) #"CH3","CH2"
     k1,k2 = _parse_kv(kk,',')
     return (string(_parse_group_string_key(k1,String)),string(_parse_group_string_key(k2,String)))
end

gc_get_comp(x::AbstractString) = x
gc_get_comp(x) = first(x)
gc_get_group(x::AbstractString) = nothing
gc_get_group(x) = last(x)


function GroupParam(gccomponents::Vector,
    group_locations::Vector{String}=String[];
    group_userlocations::Vector{String}=String[],
    verbose::Bool = false,
    grouptype = :unknown)
    options = ParamOptions(;group_userlocations,verbose)
    return GroupParam(gccomponents,group_locations,options,grouptype)
end

function GroupParam(gccomponents,
    grouplocations::Array{String,1}=String[],
    options::ParamOptions = DefaultOptions,
    grouptype = :unknown)

    # The format for gccomponents is an arary of either the species name (if it
    # available ∈ the Clapeyron database, or a tuple consisting of the species
    # name, followed by a list of group => multiplicity pairs.  For example:
    # gccomponents = ["ethane",
    #                ("hexane", ["CH3" => 2, "CH2" => 4]),
    #                ("octane", ["CH3" => 2, "CH2" => 6])]
    components = Vector{String}(undef,length(gccomponents))

    to_lookup = fill(false,length(components))
    found_gcpairs = Vector{Vector{Pair{String,Int64}}}(undef,length(gccomponents))
    components = gc_get_comp.(gccomponents)
    found_gcpairs = gc_get_group.(gccomponents)
    to_lookup = isnothing.(found_gcpairs)
    usergrouplocations = options.group_userlocations
    componentstolookup = components[to_lookup]
    filepaths = flattenfilepaths(grouplocations,usergrouplocations)

    #using parsing machinery
    if any(to_lookup)
        allparams,allnotfoundparams = createparams(componentstolookup, filepaths, options, :group) #merge all found params
        raw_result, _ = compile_params(componentstolookup,allparams,allnotfoundparams,options) #generate ClapeyronParams
        raw_groups = raw_result["groups"] #SingleParam{String}
        is_valid_param(raw_groups,options) #this will check if we actually found all params, via single missing detection.
        groupsourcecsvs = raw_groups.sourcecsvs
        
        if haskey(allparams,"groups")
            _grouptype = allparams["groups"].grouptype
        else
            _grouptype = grouptype
        end
    else
        
        _grouptype = fast_parse_grouptype(filepaths)
        if _grouptype != grouptype
            _grouptype = grouptype
        end
        groupsourcecsvs = filepaths
    end
    gccomponents_parsed = PARSED_GROUP_VECTOR_TYPE(undef,length(gccomponents))
    j = 0
    for (i,needs_to_parse_group_i) ∈ pairs(to_lookup)
        if needs_to_parse_group_i #we looked up this component, and if we are here, it exists.
            j += 1
            gcdata = _parse_group_string(raw_groups.values[j])
            gccomponents_parsed[i] = (components[i],gcdata)
        else
            gccomponents_parsed[i] = (components[i],found_gcpairs[i])
        end
    end
    return GroupParam(gccomponents_parsed,_grouptype,groupsourcecsvs)
end

function SecondOrderGroupParam(gccomponents,intragccomponents,
    grouplocations::Array{String,1}=String[],
    options::ParamOptions = DefaultOptions,
    grouptype = :unknown)


    group1 = GroupParam(gccomponents,grouplocations,options,grouptype)
    if length(intragccomponents) == 0 && length(gccomponents) != 0
        components = group1.components
        found_gcpairs = fill(nothing,length(components))
        to_lookup = fill(true,length(components))
    elseif length(intragccomponents) == length(gccomponents)
        components = gc_get_comp.(intragccomponents)
        found_gcpairs = gc_get_group.(intragccomponents)
        to_lookup = isnothing.(intragccomponents)
    else
        throw(error("GC components and intra-GC components should have the same length when provided."))
    end

    usergrouplocations = options.group_userlocations
    componentstolookup = components[to_lookup]
    filepaths = flattenfilepaths(grouplocations,usergrouplocations)

    if any(to_lookup)
        allparams,allnotfoundparams = createparams(componentstolookup, filepaths, options, :intragroup)
        raw_result, _ = compile_params(componentstolookup,allparams,allnotfoundparams,options) #generate ClapeyronParams
        raw_groups = raw_result["intragroups"] #SingleParam{String}
        is_valid_param(raw_groups,options) #this will check if we actually found all params, via single missing detection.        
        groupsourcecsvs = raw_groups.sourcecsvs
        if haskey(allparams,"intragroups")
            _grouptype = allparams["intragroups"].grouptype
        else
            _grouptype = grouptype
        end
    else
        _grouptype = fast_parse_grouptype(filepaths)
        if _grouptype != grouptype
            _grouptype = grouptype
        end
    end

    if _grouptype != group1.grouptype && _grouptype != :unkwown
        error_different_grouptype(group1.grouptype,group1.grouptype)
    end
    
    gccomponents_parsed = Vector{Pair{String,Vector{Pair{NTuple{2,String},Int}}}}(undef,length(gccomponents))
    j = 0
    for (i,needs_to_parse_group_i) ∈ pairs(to_lookup)
        if needs_to_parse_group_i #we looked up this component, and if we are here, it exists.
            j += 1
            gcdata = _parse_group_string(raw_groups.values[j],NTuple{2,String})
            gccomponents_parsed[i] = (components[i],gcdata)
        else
            gccomponents_parsed[i] = (components[i],found_gcpairs[i])
        end
    end
    return SecondOrderGroupParam(group1,gccomponents_parsed,filepaths)
end
