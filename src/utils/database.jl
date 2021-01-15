@enum CSVType singledata pairdata assocdata groupdata

"""
    getfileextension(filepoth)

A quick helper to get the file extension of any given path (without the dot).

# Examples
```julia-repl
julia> getfileextension("~/Desktop/text.txt")
"txt"
```
"""
function getfileextension(filepath::String)
    dotpos = findlast(isequal('.'), filepath)
    isnothing(dotpos) && return ""
    return filepath[dotpos+1:end]
end

"""
    getpaths(location; relativetodatabase=false)

Returns database paths that is optionally relative to OpenSAFT.jl directory.
If path is a file, then return an Array containing a single path to that file.
If path is a directory, then return an Array containing paths to all csv files in that directory.

# Examples
```julia-repl
julia> getpaths("SAFT/PCSAFT"; relativetodatabase=true)
3-element Array{String,1}:
 "/home/user/.julia/packages/OpenSAFT.jl/xxxxx/database/SAFT/PCSAFT/data_PCSAFT_assoc.csv"
 "/home/user/.julia/packages/OpenSAFT.jl/xxxxx/database/SAFT/PCSAFT/data_PCSAFT_like.csv"
 "/home/user/.julia/packages/OpenSAFT.jl/xxxxx/database/SAFT/PCSAFT/data_PCSAFT_unlike.csv"

```
"""
function getpaths(location::String; relativetodatabase::Bool=false)
    # We do not use realpath here directly because we want to make the .csv suffix optional.
    filepath = relativetodatabase ? normpath(dirname(pathof(OpenSAFT)), "..", "database", location) : location
    isfile(filepath) && return [realpath(filepath)]
    isfile(filepath * ".csv") && return [realpath(filepath * ".csv")]
    if !isdir(filepath)
        relativetodatabase ? error("The path ", location, " does not exist in the OpenSAFT database.") :
            error("The path ", location, " does not exist.")
    end
    files = joinpath.(filepath, readdir(filepath))
    return realpath.(files[isfile.(files) .& (getfileextension.(files) .== "csv")])
end

export getparams
function getparams(components::Array{String,1}, locations::Array{String,1}=String[]; userlocations::Array{String,1}=String[], asymmetricparams::Array{String,1}=String[], ignore_missingsingleparams::Bool=false, verbose::Bool=false)
    # Gets all parameters from database.
    # locations is a list of paths relative to the OpenSAFT database directory.
    # userlocations is a list of paths input by the user.
    # If parameters exist in multiple files, OpenSAFT gives priority to files in later paths.
    # asymmetricparams is a list of parameters for which matrix reflection is disabled.
    # ignore_missingsingleparams gives users the option to disable component existence check in single params.
    filepaths = string.(vcat([(getpaths.(locations; relativetodatabase=true)...)...], [(getpaths.(userlocations)...)...]))
    allcomponentsites = findsitesincsvs(components, filepaths)
    allparams, paramsourcecsvs, paramsources = createparamarrays(components, filepaths, allcomponentsites; verbose=verbose)
    return packageparams(allparams, components, allcomponentsites, paramsourcecsvs, paramsources; asymmetricparams=asymmetricparams, ignore_missingsingleparams=ignore_missingsingleparams)
end

function getparams(groups::GCParam, locations::Array{String,1}=String[]; userlocations::Array{String,1}=String[], asymmetricparams::Array{String,1}=String[], ignore_missingsingleparams::Bool=false, verbose::Bool=false)
    # For GC.
    return getparams(groups.flattenedgroups, locations; userlocations=userlocations, asymmetricparams=asymmetricparams, ignore_missingsingleparams=ignore_missingsingleparams, verbose=verbose)
end

function packageparams(allparams::Dict, components::Array{String,1}, allcomponentsites::Array{Array{String,1},1}, paramsourcecsvs::Dict{String,Set{String}}, paramsources::Dict{String,Set{String}}; asymmetricparams::Array{String,1}=String[], ignore_missingsingleparams::Bool=false)
    # Package params into their respective Structs.
    output = Dict{String,OpenSAFTParam}()
    for (param, value) ∈ allparams
        if value isa Array{<:Any,1}
            newvalue, ismissingvalues = defaultmissing(value)
            if !ignore_missingsingleparams && any(ismissingvalues)
                error("Missing values exist in single parameter ", param, ": ", value, ".")
            end
            output[param] = SingleParam(param, newvalue, ismissingvalues, components, allcomponentsites, collect(paramsourcecsvs[param]), collect(paramsources[param]))
        elseif value isa Array{<:Array,2}
            param ∉ asymmetricparams && mirrormatrix!(value)
            newvalue_ismissingvalues = defaultmissing.(value)
            newvalue = getindex.(newvalue_ismissingvalues, 1)
            ismissingvalues = getindex.(newvalue_ismissingvalues, 2)
            output[param] = AssocParam(param, newvalue, ismissingvalues, components, allcomponentsites, collect(paramsourcecsvs[param]), collect(paramsources[param]))
        elseif value isa Array{<:Any, 2}
            param ∉ asymmetricparams && mirrormatrix!(value)
            newvalue, ismissingvalues = defaultmissing(value)
            if (!ignore_missingsingleparams 
                && !all([ismissingvalues[x,x] for x ∈ 1:size(ismissingvalues,1)])
                && any([ismissingvalues[x,x] for x ∈ 1:size(ismissingvalues,1)]))
                error("Partial missing values exist in diagonal of pair parameter ", param, ": ", [value[x,x] for x ∈ 1:size(ismissingvalues,1)], ".")
            end
            output[param] = PairParam(param, newvalue, ismissingvalues, components, allcomponentsites, collect(paramsourcecsvs[param]), collect(paramsources[param]))
        else
            error("Format for ", param, " is incorrect.")
        end
    end
    return output
end

function createparamarrays(components::Array{String,1}, filepaths::Array{String,1}, allcomponentsites::Array{Array{String,1},1}; verbose::Bool=false)
    # Returns Dict with all parameters in their respective arrays.
    checkfor_clashingheaders(filepaths)
    allparams = Dict{String,Any}()
    paramsourcecsvs = Dict{String,Set{String}}()
    paramsources = Dict{String,Set{String}}()
    # Read the filepaths in reverse in order to ensure that unused sources do not get added.
    for filepath ∈ reverse(filepaths)
        csvtype = readcsvtype(filepath)
        if csvtype == groupdata
            verbose && println("Skipping groupdata csv ", filepath)
            continue
        end
        headerparams = readcsvheader(filepath)
        verbose && println("Searching for ", string(csvtype), " headers ", headerparams, " for components ", components, " at ", filepath, "...")
        foundparams, paramtypes, sources = findparamsincsv(components, filepath, headerparams; verbose=verbose)
        foundcomponents = collect(keys(foundparams))
        foundparams = swapdictorder(foundparams)
        for headerparam ∈ headerparams
            if !haskey(allparams, headerparam)
                if ismissing(paramtypes[headerparam])
                    allparams[headerparam] = createemptyparamsarray(Any, csvtype, components, allcomponentsites)
                else
                    allparams[headerparam] = createemptyparamsarray(paramtypes[headerparam], csvtype, components, allcomponentsites)
                end
            end
            if !haskey(paramsourcecsvs, headerparam)
                paramsourcecsvs[headerparam] = Set{String}()
                paramsources[headerparam] = Set{String}()
            end
            if csvtype == singledata
                isempty(foundparams) && continue
                for (component, value) ∈ foundparams[headerparam]
                    currenttype = nonmissingtype(eltype(allparams[headerparam]))
                    if paramtypes[headerparam] != currenttype
                        if currenttype <: Number && paramtypes[headerparam] <: String || currenttype <: String && paramtypes[headerparam] <: Number 
                            allparams[headerparam] = convert(Array{Any}, allparams[headerparam])
                        else
                            try
                                allparams[headerparam] = convert(Array{Union{Missing,paramtypes[headerparam]}}, allparams[headerparam])
                            catch e
                            end
                        end
                    end
                    idx = findfirst(isequal(component), components)
                    if allparams[headerparam] isa Array{<:Any,2}
                        !ismissing(allparams[headerparam][idx,idx]) && continue
                        allparams[headerparam][idx,idx] = value
                    else
                        !ismissing(allparams[headerparam][idx]) && continue
                        allparams[headerparam][idx] = value
                    end
                    push!(paramsourcecsvs[headerparam], filepath)
                    !ismissing(sources[component]) && push!(paramsources[headerparam], sources[component])
                end
            end
            if csvtype == pairdata
                if allparams[headerparam] isa Array{<:Any,1}
                    allparams[headerparam] = convertsingletopair(allparams[headerparam], outputmissing=true)
                end
                isempty(foundparams) && continue
                for (componentpair, value) ∈ foundparams[headerparam]
                    currenttype = nonmissingtype(eltype(allparams[headerparam]))
                    if paramtypes[headerparam] != currenttype
                        if currenttype <: Number && paramtypes[headerparam] <: String || currenttype <: String && paramtypes[headerparam] <: Number 
                            allparams[headerparam] = convert(Array{Any}, allparams[headerparam])
                        else
                            try
                                allparams[headerparam] = convert(Array{Union{Missing,paramtypes[headerparam]}}, allparams[headerparam])
                            catch e 
                            end
                        end
                    end
                    idx1 = findfirst(isequal(componentpair[1]), components)
                    idx2 = findfirst(isequal(componentpair[2]), components)
                    !ismissing(allparams[headerparam][idx1,idx2]) && continue
                    allparams[headerparam][idx1,idx2] = value
                    push!(paramsourcecsvs[headerparam], filepath)
                    !ismissing(sources[componentpair]) && push!(paramsources[headerparam], sources[componentpair])
                end
            end
            if csvtype == assocdata
                isempty(foundparams) && continue
                for (assocpair, value) ∈ foundparams[headerparam]
                    currenttype = nonmissingtype(eltype(first(allparams[headerparam])))
                    if currenttype != paramtypes[headerparam]
                        if currenttype <: Number && paramtypes[headerparam] <: String || currenttype <: String && paramtypes[headerparam] <: Number 
                            allparams[headerparam] = convert(Array{Array{Any,2},2}, allparams[headerparam])
                        else
                            try
                                allparams[headerparam] = convert(Array{Array{Union{Missing,paramtypes[headerparam]},2},2}, allparams[headerparam])
                            catch e
                            end
                        end
                    end
                    idx1 = findfirst(isequal(assocpair[1][1]), components)
                    idx2 = findfirst(isequal(assocpair[1][2]), components)
                    idx21 = findfirst(isequal(assocpair[2][1]), allcomponentsites[idx1])
                    idx22 = findfirst(isequal(assocpair[2][2]), allcomponentsites[idx2])
                    !ismissing(allparams[headerparam][idx1,idx2][idx21,idx22]) && continue
                    allparams[headerparam][idx1,idx2][idx21,idx22] = value
                    push!(paramsourcecsvs[headerparam], filepath)
                    !ismissing(sources[assocpair]) && push!(paramsources[headerparam], sources[assocpair])
                end
            end
        end
    end
    return allparams, paramsourcecsvs, paramsources
end

function defaultmissing(array::Array; defaultvalue=nothing)
    # Returns a non-missing Array with missing values replaced by default values.
    # Also returns an Array of Bools to retain overwritten information.
    arraycopy = deepcopy(array)
    type = nonmissingtype(eltype(array))
    ismissingvalues = ismissing.(array)
    if isnothing(defaultvalue)
        if type == Union{} #  If it only contains missing
            type = Any
        elseif type == Any
            # This means that it has a mix of Number and String types.
            # If that was the case, convert it all to String.
            arraycopy[ismissingvalues] .= ""
            arraycopy = string.(arraycopy)
            type = String
        elseif type <: Number
            arraycopy[ismissingvalues] .= 0
        elseif type <: AbstractString
            arraycopy[ismissingvalues] .= ""
        else
            error("Unsupported type.")
        end
    else
        arraycopy[ismissingvalues] .= defaultvalue
    end
    return convert(Array{type}, arraycopy), convert(Array{Bool}, ismissingvalues)
end

function swapdictorder(dict::Dict)
    # Swap the first two levels in a nested dictionary.
    isempty(dict) && return dict
    output = Dict()
    outerkeys = keys(dict)
    try
        innerkeys = keys(first(values(dict)))
        for innerkey ∈ innerkeys, outerkey ∈ outerkeys
            if !haskey(output, innerkey)
                output[innerkey] = Dict{Any,Any}(outerkey => dict[outerkey][innerkey])
            end
            push!(output[innerkey], outerkey => dict[outerkey][innerkey])
        end
    catch e
        error("The format of the nested dictionary is not correct.")
    end
    return output
end

function findparamsincsv(components::Array{String,1}, filepath::String, headerparams::Array{String,1}; columnreference::String="species", sitecolumnreference::String="site", sourcecolumnreference::String="source", verbose::Bool=false, ignore_missingsingleparams::Bool=false)
    # Returns a Dict with all matches in a particular file for one parameter.
    normalised_columnreference = normalisestring(columnreference)
    csvtype = readcsvtype(filepath)
    df = CSV.File(filepath; header=3, pool=0)
    csvheaders = String.(Tables.columnnames(df))
    normalised_csvheaders = normalisestring.(csvheaders)
    normalised_headerparams = normalisestring.(headerparams)
    normalised_headerparams ⊈ normalised_csvheaders && error("Headers ", setdiff(normalised_headerparams, normalised_csvheaders), " not present in csv header.")

    foundvalues = Dict()
    paramtypes = Dict(headerparams .=> [Tables.columntype(df, Symbol(x)) for x ∈ headerparams])
    sources = Dict()

    normalised_sourcecolumnreference = normalisestring(sourcecolumnreference)
    getsources = false
    if normalised_sourcecolumnreference ∈ normalised_csvheaders
        getsources = true
        sourcecolumn = Symbol(normalised_csvheaders[findfirst(isequal(normalised_sourcecolumnreference), normalised_csvheaders)])
    end

    if csvtype == singledata
        lookupcolumnindex = findfirst(isequal(normalised_columnreference), normalised_csvheaders)
        isnothing(lookupcolumnindex) && error("Header ", normalised_columnreference, " not found.")
        lookupcolumn = Symbol(csvheaders[lookupcolumnindex])
        for row ∈ Tables.rows(df)
            component = row[lookupcolumn]
            component ∉ components && continue
            verbose && print("Found component: ", component)
            foundvalues[component] = Dict()
            for headerparam ∈ headerparams
                foundvalues[component][headerparam] = row[Symbol(headerparam)]
            end
            verbose && println(" with values ", foundvalues[component])
            if getsources
                source = row[sourcecolumn]
                sources[component] = source
            else
                sources[component] = missing
            end
        end
    elseif csvtype == pairdata
        lookupcolumnindex1 = findfirst(x -> x[1:end-1] == normalised_columnreference && x[end] == '1', normalised_csvheaders)
        lookupcolumnindex2 = findfirst(x -> x[1:end-1] == normalised_columnreference && x[end] == '2', normalised_csvheaders)
        isnothing(lookupcolumnindex1) && error("Header ", normalised_columnreference * "1", " not found.")
        isnothing(lookupcolumnindex2) && error("Header ", normalised_columnreference * "2", " not found.")
        lookupcolumn1 = Symbol(csvheaders[lookupcolumnindex1])
        lookupcolumn2 = Symbol(csvheaders[lookupcolumnindex2])
        for row ∈ Tables.rows(df)
            component1 = row[lookupcolumn1]
            component2 = row[lookupcolumn2]
            (component1 ∉ components || component2 ∉ components) && continue
            componentpair = (component1, component2)
            verbose && print("Found component pair: ", componentpair)
            foundvalues[componentpair] = Dict()
            for headerparam ∈ headerparams
                foundvalues[componentpair][headerparam] = row[Symbol(headerparam)]
            end
            verbose && println(" with values ", foundvalues[componentpair])
            if getsources
                source = row[sourcecolumn]
                sources[componentpair] = source
            else
                sources[componentpair] = missing
            end
        end
    elseif csvtype == assocdata
        lookupcolumnindex1 = findfirst(x -> x[1:end-1] == normalised_columnreference && x[end] == '1', normalised_csvheaders)
        lookupcolumnindex2 = findfirst(x -> x[1:end-1] == normalised_columnreference && x[end] == '2', normalised_csvheaders)
        isnothing(lookupcolumnindex1) && error("Header ", normalised_columnreference * "1", " not found.")
        isnothing(lookupcolumnindex2) && error("Header ", normalised_columnreference * "2", " not found.")
        lookupcolumn1 = Symbol(csvheaders[lookupcolumnindex1])
        lookupcolumn2 = Symbol(csvheaders[lookupcolumnindex2])

        normalised_sitecolumnreference = normalisestring(sitecolumnreference)
        lookupsitecolumnindex1 = findfirst(x -> x[1:end-1] == normalised_sitecolumnreference && x[end] == '1', normalised_csvheaders)
        lookupsitecolumnindex2 = findfirst(x -> x[1:end-1] == normalised_sitecolumnreference && x[end] == '2', normalised_csvheaders)
        isnothing(lookupcolumnindex1) && error("Header ", normalised_sitecolumnreference * "1", " not found.")
        isnothing(lookupcolumnindex2) && error("Header ", normalised_sitecolumnreference * "2", " not found.")
        lookupsitecolumn1 = Symbol(csvheaders[lookupsitecolumnindex1])
        lookupsitecolumn2 = Symbol(csvheaders[lookupsitecolumnindex2])

        for row ∈ Tables.rows(df)
            component1 = row[lookupcolumn1]
            component2 = row[lookupcolumn2]
            (component1 ∉ components || component2 ∉ components) && continue
            site1 = row[lookupsitecolumn1]
            site2 = row[lookupsitecolumn2]
            assocpair = ((component1, component2), (site1, site2))
            verbose && print("Found assoc pair: ", assocpair)
            foundvalues[assocpair] = Dict()
            for headerparam ∈ headerparams
                foundvalues[assocpair][headerparam] = row[Symbol(headerparam)]
            end
            verbose && println(" with values ", foundvalues[assocpair])
            if getsources
                source = row[sourcecolumn]
                sources[assocpair] = source
            else
                sources[assocpair] = missing
            end
        end
    else
        error("File is of type ", String(csvtype), " and cannot be read with this function.")
    end
    return foundvalues, paramtypes, sources
end

function normalisestring(str::String)
    return lowercase(replace(str, r"[ \-\_]" => ""))
end

function findgroupsincsv(components::Array{String,1}, filepath::String; columnreference::String="species", groupcolumnreference::String="groups", verbose::Bool=false)
    # Returns a Dict with the group string that will be parsed in buildspecies.
    csvtype = readcsvtype(filepath)
    csvtype != groupdata && return Dict{String,String}()
    normalised_columnreference = normalisestring(columnreference)
    normalised_groupcolumnreference = normalisestring(groupcolumnreference)
    csvtype = readcsvtype(filepath)
    df = CSV.File(filepath; header=3)
    csvheaders = String.(Tables.columnnames(df))
    normalised_csvheaders = normalisestring.(csvheaders)

    normalised_columnreference ∉ normalised_csvheaders && error("Header ", normalised_columnreference, " not found.")
    normalised_groupcolumnreference ∉ normalised_csvheaders && error("Header ", normalised_groupcolumnreference, " not found.")

    foundgroups = Dict{String,String}()

    lookupcolumnindex = findfirst(isequal(normalised_columnreference), normalised_csvheaders)
    isnothing(lookupcolumnindex) && error("Header ", normalised_columnreference, " not found.")
    lookupcolumn = Symbol(csvheaders[lookupcolumnindex])

    lookupgroupcolumnindex = findfirst(isequal(normalised_groupcolumnreference), normalised_csvheaders)
    isnothing(lookupgroupcolumnindex) && error("Header ", normalised_groupcolumnreference, " not found.")
    lookupgroupcolumn = Symbol(csvheaders[lookupgroupcolumnindex])
    for row ∈ Tables.rows(df)
        component = row[lookupcolumn]
        component ∉ components && continue
        verbose && print("Found component: ", component)
        foundgroups[component] = row[lookupgroupcolumn]
        verbose && println(" with groups ", foundgroups[component])
    end
    return foundgroups
end

function readcsvtype(filepath::String)
    # Searches for type from second line of CSV.
    keywords = ["like", "single", "unlike", "pair", "assoc", "group", "groups"]
    words = split(lowercase(rstrip(getline(filepath, 2), ',')), ' ')
    foundkeywords = intersect(words, keywords)
    isempty(foundkeywords) && error("Unable to determine type of database", filepath, ". Check that keyword is present on Line 2.")
    length(foundkeywords) > 1 && error("Multiple keywords found in database ", filepath, ": ", foundkeywords)
    first(foundkeywords) == "single" && return singledata
    first(foundkeywords) == "like" && return singledata
    first(foundkeywords) == "pair" && return pairdata
    first(foundkeywords) == "unlike" && return pairdata
    first(foundkeywords) == "assoc" && return assocdata
    first(foundkeywords) == "group" && return groupdata
    first(foundkeywords) == "groups" && return groupdata
end

function getline(filepath::String, selectedline::Int)
    # Simple function to return text from filepath at selectedline.
    open(filepath) do file
        linecount = 1
        for line ∈ eachline(file)
            linecount == selectedline && return line
            linecount += 1
        end
        error("Selected line number exceeds number of lines in file")
    end
end
            
function readcsvheader(filepath::String; headerline::Int = 3)
    # Returns array of filtered header strings at line 3.
    headers = split(getline(filepath, headerline), ',')
    ignorelist = ["source", "species", "dipprnumber", "smiles", "site"]
    return String.(filter(x -> replace.(lowercase(x), r"[ \-\_\d]" => "") ∉ ignorelist, headers))
end

function checkfor_clashingheaders(filepaths::Array{String,1})
    # Raises an error if the header of any assoc parameter clashes with a non-assoc parameter.
    headerparams = []
    headerparams_assoc = []
    for filepath in filepaths
        csvtype = readcsvtype(filepath)
        if csvtype == singledata || csvtype == pairdata
            append!(headerparams, readcsvheader(filepath))
        elseif csvtype == assocdata
            append!(headerparams_assoc, readcsvheader(filepath))
        end
    end
    clashingheaders = intersect(headerparams, headerparams_assoc)
    !isempty(clashingheaders) && error("Headers ", clashingheaders, " appear in both loaded assoc and non-assoc files.")
end

function findsitesincsvs(components::Array{String,1}, filepaths::Array{String,1}; columnreference::String="species", sitecolumnreference::String="site", verbose::Bool=false)
    # Look for all relevant sites in the database.
    # Note that this might not necessarily include all sites associated with a component.
    normalised_columnreference = normalisestring(columnreference)
    normalised_sitecolumnreference = normalisestring(sitecolumnreference)
    sites = Dict(components .=> [Set{String}() for _ ∈ 1:length(components)])
    for filepath ∈ filepaths
        csvtype = readcsvtype(filepath)
        csvtype != assocdata && continue

        df = CSV.File(filepath; header=3)
        csvheaders = String.(Tables.columnnames(df))
        normalised_csvheaders = normalisestring.(String.(Tables.columnnames(df)))

        lookupcolumnindex1 = findfirst(x -> x[1:end-1] == normalised_columnreference && x[end] == '1', normalised_csvheaders)
        lookupcolumnindex2 = findfirst(x -> x[1:end-1] == normalised_columnreference && x[end] == '2', normalised_csvheaders)
        isnothing(lookupcolumnindex1) && error("Header ", normalised_columnreference * "1", " not found.")
        isnothing(lookupcolumnindex2) && error("Header ", normalised_columnreference * "2", " not found.")
        lookupcolumn1 = Symbol(csvheaders[lookupcolumnindex1])
        lookupcolumn2 = Symbol(csvheaders[lookupcolumnindex2])

        normalised_sitecolumnreference = normalisestring(sitecolumnreference)
        lookupsitecolumnindex1 = findfirst(x -> x[1:end-1] == normalised_sitecolumnreference && x[end] == '1', normalised_csvheaders)
        lookupsitecolumnindex2 = findfirst(x -> x[1:end-1] == normalised_sitecolumnreference && x[end] == '2', normalised_csvheaders)
        isnothing(lookupcolumnindex1) && error("Header ", normalised_sitecolumnreference * "1", " not found.")
        isnothing(lookupcolumnindex2) && error("Header ", normalised_sitecolumnreference * "2", " not found.")
        lookupsitecolumn1 = Symbol(csvheaders[lookupsitecolumnindex1])
        lookupsitecolumn2 = Symbol(csvheaders[lookupsitecolumnindex2])

        for row ∈ Tables.rows(df)
            component1 = row[Symbol(lookupcolumn1)]
            component2 = row[Symbol(lookupcolumn2)]
            if component1 ∈ components && component2 ∈ components
                push!(sites[component1], row[Symbol(lookupsitecolumn1)])
                push!(sites[component2], row[Symbol(lookupsitecolumn2)])
            end
        end
    end
    output = Array{Array{String,1}}(undef, 0)
    for component ∈ components
        push!(output, collect(sites[component]))
    end
    verbose && println("Found sites for ", components, " are ", output, ".")
    return output
end

function createemptyparamsarray(datatype::Type, csvtype::CSVType, components::Array{String,1})
    allcomponentsites = [String[] for _ in 1:length(components)]
    return createemptyparamsarray(datatype, csvtype, components, allcomponentsites)
end

function createemptyparamsarray(datatype::Type, csvtype::CSVType, components::Array{String,1}, allcomponentsites::Array{Array{String,1},1})
    # Creates a missing array of the appropriate size.
    componentslength = length(components)
    csvtype == singledata && return (Array{Union{Missing,datatype}}(undef, componentslength) .= missing)
    csvtype == pairdata && return (Array{Union{Missing,datatype}}(undef, componentslength, componentslength) .= missing)
    if csvtype == assocdata
        output = Array{Array{Union{Missing,datatype},2},2}(undef, componentslength, componentslength)
        for i ∈ 1:componentslength, j ∈ 1:componentslength
            output[i,j] = (Array{Union{Missing,datatype}}(undef, length(allcomponentsites[i]), length(allcomponentsites[j])) .= missing)
        end
        return output
    end
end

function createemptyparamsarray(csvtype::CSVType, components::Array{String,1})
    allcomponentsites = [String[] for _ in 1:length(components)]
    return createemptyparamsarray(Any, csvtype, components, allcomponentsites)
end

function createemptyparamsarray(csvtype::CSVType, components::Array{String,1}, allcomponentsites::Array{Array{String,1},1})
    return createemptyparamsarray(Any, csvtype, components, allcomponentsites)
end

function convertsingletopair(params::Array{T,1}; outputmissing::Bool=false) where T
    # Returns a missing square matrix with its diagonal matrix replaced by the given parameters. 
    paramslength = length(params)
    if outputmissing
        output = Array{Union{Missing,T}}(undef, paramslength, paramslength) .= missing
    else
        if T == Any
            output = Array{T}(undef, paramslength, paramslength) .= 0
        elseif T <: Number
            output = Array{T}(undef, paramslength, paramslength) .= 0
        elseif T <: AbstractString
            output = Array{T}(undef, paramslength, paramslength) .= ""
        else
            error("Data type ", T, " not supported.")
        end
    end
    for i = 1:paramslength
        output[i,i] = params[i]
    end
    return output
end

function mirrormatrix!(matrix::Array{T,2}) where T
    # Mirrors a square matrix.
    matrixsize = size(matrix)
    matrixsize[1] != matrixsize[2] && error("Matrix is not square.")
    for i ∈ 2:matrixsize[1], j ∈ 1:i-1
        lowervalue = matrix[i,j]
        uppervalue = matrix[j,i]
        if !ismissing(lowervalue) && !ismissing(uppervalue) && lowervalue != uppervalue
            error("Dissimilar non-zero entries exist across diagonal:", lowervalue, ", ", uppervalue)
        end
        !ismissing(lowervalue) && (matrix[j,i] = lowervalue)
        !ismissing(uppervalue) && (matrix[i,j] = uppervalue)
    end
    return matrix
end

function mirrormatrix!(matrix::Array{Array{T,2},2}) where T
    # Mirrors a square matrix.
    matrixsize = size(matrix)
    matrixsize[1] != matrixsize[2] && error("Matrix is not square.")
    [mirrormatrix!(matrix[i,i]) for i ∈ 1:matrixsize[1]]
    for i ∈ 2:matrixsize[1], j ∈ 1:i-1
        matrixsize2 = size(matrix[i,j])
        for a ∈ 1:matrixsize2[1], b ∈ 1:matrixsize2[2]
            lowervalue = matrix[i,j][a,b]
            uppervalue = matrix[j,i][b,a]
            if !ismissing(lowervalue) && !ismissing(uppervalue) && lowervalue != uppervalue
                error("Dissimilar non-zero entries exist across diagonal:", lowervalue, ", ", uppervalue)
            end
            !ismissing(lowervalue) && (matrix[j,i][b,a] = lowervalue)
            !ismissing(uppervalue) && (matrix[i,j][a,b] = uppervalue)
        end
    end
    return matrix
end


function buildspecies(gccomponents::Array{<:Any,1}, grouplocations::Array{String,1}=String[]; usergrouplocations::Array{String,1}=String[], verbose::Bool=false)
    # The format for gccomponents is an arary of either the species name (if it
    # available in the OpenSAFT database, or a tuple consisting of the species
    # name, followed by a list of group => multiplicity pairs.  For example:
    # gccomponents = ["ethane",
    #                ("hexane", ["CH3" => 2, "CH2" => 4]),
    #                ("octane", ["CH3" => 2, "CH2" => 6])]
    BuildSpeciesType = Union{Tuple{String, Array{Pair{String, Int64},1}}, String, Tuple{String}}
    any(.!(isa.(gccomponents, BuildSpeciesType))) && error("The format of the components is incorrect.")
    filepaths = string.(vcat([(getpaths.(grouplocations; relativetodatabase=true)...)...], [(getpaths.(usergrouplocations)...)...]))
    components = String[]
    allcomponentgroups = Array{Array{String,1},1}(undef, 0)
    allcomponentngroups = Array{Array{Int,1},1}(undef, 0)
    componentstolookup = String[]
    append!(componentstolookup, [x for x ∈ gccomponents[isa.(gccomponents, String)]])
    append!(componentstolookup, [first(x) for x ∈ gccomponents[isa.(gccomponents, Tuple{String})]])
    allfoundcomponentgroups = Dict{String,String}()
    groupsourcecsvs = String[]
    for filepath in filepaths
        csvtype = readcsvtype(filepath)
        if csvtype != groupdata
            verbose && println("Skipping ", csvtype, " csv at ", filepath)
            continue
        end
        verbose && println("Searching for groups for components ", componentstolookup, " at ", filepath, "...")
        merge!(allfoundcomponentgroups, findgroupsincsv(componentstolookup, filepath; verbose=verbose))
        append!(groupsourcecsvs, [filepath])
    end
    for gccomponent ∈ gccomponents
        if gccomponent isa Tuple{String, Array{Pair{String, Int64},1}}
            append!(components, [gccomponent[1]])
            groupsandngroups = gccomponent[2]
        elseif gccomponent isa String
            !haskey(allfoundcomponentgroups, gccomponent) && error("Predefined component ", gccomponent, " not found in any group input csvs.")
            append!(components, [gccomponent])
            groupsandngroups = eval(Meta.parse(allfoundcomponentgroups[gccomponent]))
        elseif gccomponent isa Union{String, Tuple{String}}
            !haskey(allfoundcomponentgroups, first(gccomponent)) && error("Predefined component ", gccomponent, " not found in any group input csvs.")
            append!(components, [first(gccomponent)])
            groupsandngroups = eval(Meta.parse(allfoundcomponentgroups[first(gccomponent)]))
        end
        componentgroups = String[]
        componentngroups = Int[]
        for (componentgroup, ncomponentgroup) ∈ groupsandngroups
            append!(componentgroups, [componentgroup])
            append!(componentngroups, [ncomponentgroup])
        end
        append!(allcomponentgroups, [componentgroups])
        append!(allcomponentngroups, [componentngroups])
    end
    flattenedgroups = unique([(allcomponentgroups...)...])
    allcomponentnflattenedgroups = [zeros(Int, length(flattenedgroups)) for _ in 1:length(components)]
    for i in 1:length(components)
        for (k, group) in enumerate(flattenedgroups)
            if group in allcomponentgroups[i]
                allcomponentnflattenedgroups[i][k] = allcomponentngroups[i][findfirst(isequal(group), allcomponentgroups[i])]
            else
                allcomponentnflattenedgroups[i][k] = 0
            end
        end
    end
    return GCParam(components, allcomponentgroups, allcomponentngroups, flattenedgroups, allcomponentnflattenedgroups, groupsourcecsvs)
end
