function fast_parse_grouptype(filepaths::Vector{String})
    #only parses grouptype, if present in any CSV, is used. if not, return unknown
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
                if new_grouptype != :unknown #for backwards compatibility
                    error_different_grouptype(grouptype,new_grouptype)
                end
            end
        end
    end
    grouptype == :not_set && (grouptype = :unknown)
    return grouptype
end

#used to parse """["CH => 2, "OH" = 2, "CH3" => 2]"""
function _parse_group_string(gc::String,gctype=String)
    gc_strip = strip(gc)
    if startswith(gc_strip,"[") && endswith(gc_strip,"]")
        gc_without_brackets = chop(gc_strip,head = 1,tail = 1)
        if gctype == String
            gcpairs = split(gc_without_brackets,",")
        else
            gc_without_brackets = replace(gc_without_brackets,",(" => ";;(")
            gc_without_brackets = replace(gc_without_brackets,", (" => ";;(")
            gc_without_brackets = replace(gc_without_brackets,",  (" => ";;(")

            gcpairs = split(gc_without_brackets,";;")
        end
        result = Vector{Pair{gctype,Int}}(undef,length(gcpairs))
        for (i,gcpair) ∈ pairs(gcpairs)
            raw_group_i,raw_num = _parse_kv(gcpair,"=>")
            if (startswith(raw_group_i,"\"") && endswith(raw_group_i,"\"")) ||
                (startswith(raw_group_i,"(") && endswith(raw_group_i,")"))
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

#getting the component name part
gc_get_comp(x::AbstractString) = x
gc_get_comp(x) = first(x)

#getting the group part
gc_get_group(x::AbstractString) = nothing
gc_get_group(x::Tuple{Any,Any}) = x[2] #first index is the component, second is the group, third is the bond info.
gc_get_group(x::Tuple{Any,Any,Any}) = x[2]
gc_get_group(x::Pair) = last(x)

#getting the intragroup part
gc_get_intragroup(x::AbstractString) = nothing
gc_get_intragroup(x::Tuple{Any,Any}) = nothing
gc_get_intragroup(x::Tuple{Any,Any,Any}) = x[3]
gc_get_intragroup(x::Pair) = nothing


function GroupParam(gccomponents::Union{AbstractVector,Tuple{String,Vector},Pair{String,Vector}},
    group_locations=String[];
    group_userlocations = String[],
    verbose::Bool = false,
    grouptype = :unknown)
    options = ParamOptions(;group_userlocations=String.(group_userlocations),verbose)
    if isa(gccomponents,AbstractVector)
        return GroupParam(gccomponents,String.(group_locations),options,grouptype)
    else
        return GroupParam([gccomponents],String.(group_locations),options,grouptype)
    end
end

function GroupParam(gccomponents,
    grouplocations=String[],
    options::ParamOptions = DefaultGroupOptions,
    grouptype = :unknown)
    # The format for gccomponents is an arary of either the species name (if it
    # available ∈ the Clapeyron database, or a tuple consisting of the species
    # name, followed by a list of group => multiplicity pairs.  For example:
    # gccomponents = ["ethane",
    #                ("hexane", ["CH3" => 2, "CH2" => 4]),
    #                ("octane", ["CH3" => 2, "CH2" => 6])]

    _gccomponents = format_gccomponents(gccomponents)
    components = map(gc_get_comp,_gccomponents)
    found_gcpairs = map(gc_get_group,_gccomponents)
    found_intragcpairs = map(gc_get_intragroup,_gccomponents)
    return __GroupParam(components,found_gcpairs,found_intragcpairs,grouplocations,options,grouptype)
end

function __GroupParam(components,found_gcpairs,found_intragcpairs,grouplocations,options,grouptype)

    to_lookup = isnothing.(found_gcpairs)
    usergrouplocations = String.(options.group_userlocations)
    componentstolookup = components[to_lookup]
    filepaths = flattenfilepaths(grouplocations,usergrouplocations)

    gccomponents_parsed = PARSED_GROUP_VECTOR_TYPE(undef,length(components))

    #fill gccomponents_parsed
    for i in 1:length(components)
        if !to_lookup[i]
            gccomponents_parsed[i] = (components[i],found_gcpairs[i])
        end
    end

    #using parsing machinery
    if any(to_lookup)
        allparams,allnotfoundparams = createparams(componentstolookup, filepaths, options, :group) #merge all found params
        raw_result = compile_params(componentstolookup,allparams,allnotfoundparams,nothing,options) #generate ClapeyronParams
        raw_groups = raw_result["groups"] #SingleParam{String}

        is_valid_param(raw_groups,options) #this will check if we actually found all params, via single missing detection.
        groupsourcecsvs = raw_groups.sourcecsvs

        if haskey(allparams,"groups")
            _grouptype = allparams["groups"].grouptype
        else
            _grouptype = grouptype
        end

        j = 0
        for (i,needs_to_parse_group_i) ∈ pairs(to_lookup)
            if needs_to_parse_group_i #we looked up this component, and if we are here, it exists.
                j += 1
                gcdata = _parse_group_string(raw_groups.values[j])
                gccomponents_parsed[i] = (components[i],gcdata)
            end
        end
    else
        _grouptype = fast_parse_grouptype(filepaths)
        if _grouptype != grouptype && grouptype != :unknown
            _grouptype = grouptype
        end
        groupsourcecsvs = filepaths
        raw_result = Dict{ClapeyronParam,String}()
    end
    #group without second order information
    group1 = GroupParam(gccomponents_parsed,_grouptype,groupsourcecsvs,nothing)

    #if we don't have second order information, return early
    if all(isnothing,found_intragcpairs) && !haskey(raw_result,"intragroups")
        return group1
    end

    #parse second order information
    empty_structgc = Vector{Pair{NTuple{2,String},Int}}(undef,0)
    structgc_components_parsed = fill(empty_structgc,length(components))
    to_lookup_intragroups = isnothing.(found_intragcpairs)

    for i in 1:length(components)
        gcpair_i = found_intragcpairs[i]
        if !to_lookup_intragroups[i]
            structgc_components_parsed[i] = gcpair_i
        else
            #if the value is missing, try looking up if its connectivity is trivial
            trivial_intragroup = build_trivial_intragroup(group1,i,to_lookup_intragroups[i])
            if !isnothing(trivial_intragroup)
                structgc_components_parsed[i] = trivial_intragroup
                to_lookup_intragroups[i] = false
            end
        end
    end

    #parse groups from database if not already set
    if any(to_lookup_intragroups)
        k = 0
        intragroups = get(raw_result,"intragroups") do
            SingleParam("intragroups",components)
        end
        for i ∈ 1:length(components)
            looked_up = to_lookup[i]
            needs_to_parse_intragroup_i = to_lookup_intragroups[i]
            #if the value is not looked up, the preparsed result stays empty
            if looked_up && needs_to_parse_intragroup_i
                k = count(@view(to_lookup[1:i]))
                if !intragroups.ismissingvalues[k]
                    structgc_components_parsed[i] = _parse_group_string(intragroups.values[k],NTuple{2,String})
                end
            end

        end
    end
    #add intragroup data to the current group param, it takes care of missing values.
    build_gc_intragroups!(group1,structgc_components_parsed)
    return group1
end

function build_trivial_intragroup(groups::GroupParam,i,lookup)
    lookup || return nothing
    groupnames = groups.groups[i]
    n = groups.n_groups[i]
    len = length(groupnames)
    len > 2 && return nothing
    if len == 1
        ni = only(n)
        if ni == 1 || ni == 2
            nij = ni - 1
            return [(groupnames[1],groupnames[1]) => nij]
        else
            return nothing
        end
    else #len == 2
        if n[1] == n[2] == 1
            return [(groupnames[1],groupnames[2]) => 1]
        else
            return nothing
        end
    end
end
