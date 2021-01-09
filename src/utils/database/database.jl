using CSV, Tables
include("OpenSAFTParam.jl")
include("visualisation.jl")

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
    getdatabasepaths(model)

Returns database paths relative to OpenSAFT.jl directory.
If path is a file, then return an Array containing a single path to that file.
If path is a directory, then return an Array containing paths to all csv files in that directory.

# Examples
```julia-repl
julia> getdatabasepaths("SAFT/PCSAFT")
3-element Array{String,1}:
 "/home/user/.julia/packages/OpenSAFT.jl/xxxxx/src/../database/SAFT/PCSAFT/data_PCSAFT_assoc.csv"
 "/home/user/.julia/packages/OpenSAFT.jl/xxxxx/src/../database/SAFT/PCSAFT/data_PCSAFT_like.csv"
 "/home/user/.julia/packages/OpenSAFT.jl/xxxxx/src/../database/SAFT/PCSAFT/data_PCSAFT_unlike.csv"

```
"""
function getdatabasepaths(model::String)
    filepath = joinpath(dirname(pathof(OpenSAFT)), "../database", model)
    isfile(filepath) && return [filepath]
    isfile(filepath * ".csv") && return [filepath * ".csv"]
    !isdir(filepath) && error("The directory ", model, " does not exist in the OpenSAFT database.")
    files = joinpath.(filepath, readdir(filepath))
    return files[isfile.(files) .& (getfileextension.(files) .== "csv")]
end

"""
    getuserpaths(model)

Same as the above, but it directly takes the (relative) path that is given as input.
```
"""
function getuserpaths(model::String)
    filepath = model
    isfile(filepath) && return [filepath]
    isfile(filepath * ".csv") && return [filepath * ".csv"]
    !isdir(filepath) && error("The directory ", filepath, " does not exist.")
    files = joinpath.(filepath, readdir(filepath))
    return files[isfile.(files) .& (getfileextension.(files) .== "csv")]
end

function getmodelname(models::Array{String,1}, usermodels::Array{String,1})
    # Try to guess the name of the model.
    # It will take the name of the last given directory, checking models before usermodels.
    # It's not foolproof, so it is highly recommended that you name your model.
    if !isempty(models)
        for model ∈ reverse(models)
            filepath = joinpath(dirname(pathof(OpenSAFT)), "../database", model)
            isdir(filepath) && return basename(filepath)
        end
    end
    if !isempty(usermodels)
        for usermodel ∈ reverse(usermodels)
            filepath = usermodel
            isdir(filepath) && return basename(filepath)
        end
    end
    return "unnamed"
end

function getparams(components::Array{String,1}, models::Array{String,1}=String[]; usermodels::Array{String,1}=String[], modelname="", asymmetric_pairparams::Array{String,1}=String[], ignore_missingsingleparams=false, verbose=false)
    # Gets all parameters from database.
    # models is a list of paths relative to the OpenSAFT database directory.
    # usermodels is a list of paths input by the user.
    # If parameters exist in multiple files, OpenSAFT gives priority to files in later paths.
    # asymmetric_pairparams is a list of parameters for which matrix reflection is disabled.
    # ignore_missingsingleparams gives users the option to disable component existence check in single params.
    filepaths = string.(vcat([(getdatabasepaths.(models)...)...], [(getuserpaths.(usermodels)...)...]))
    allcomponentsites = findsitesincsvs(filepaths, components)
    allparams, paramsources = createparamarrays(filepaths, components, allcomponentsites; verbose=verbose)
    if modelname == ""
        modelname = getmodelname(models, usermodels)
    end
    finaldict = packageparams(allparams, components, allcomponentsites, paramsources, modelname; asymmetric_pairparams=asymmetric_pairparams, ignore_missingsingleparams=ignore_missingsingleparams)
    return finaldict
end

function getparams(groups::GCParam, models::Array{String,1}=String[]; usermodels::Array{String,1}=String[], modelname="", asymmetric_pairparams::Array{String,1}=String[], ignore_missingsingleparams=false, verbose=false)
    # For GC.
    return getparams(groups.flattenedgroups, models; usermodels=usermodels, modelname=modelname, asymmetric_pairparams=asymmetric_pairparams, ignore_missingsingleparams=ignore_missingsingleparams, verbose=verbose)
end

function packageparams(allparams::Dict, components::Array{String,1}, allcomponentsites::Array{Array{String,1},1}, paramsources::Dict{String,Set{String}}, modelname::String; asymmetric_pairparams::Array{String,1}=String[], ignore_missingsingleparams=false)
    # Package params into their respective Structs.
    output = Dict{String,OpenSAFTParam}()
    for (param, value) ∈ allparams
        if value isa Array{<:Any,1}
            newvalue, ismissingvalues = defaultmissing(value)
            if !ignore_missingsingleparams && any(ismissingvalues)
                error("Missing values exist in single parameter ", param, ": ", value, ".")
            end
            output[param] = SingleParam(param, newvalue, ismissingvalues, components, allcomponentsites, modelname, collect(paramsources[param]))
        elseif value isa Array{<:Array,2}
            newvalue_ismissingvalues = defaultmissing.(value)
            newvalue = getindex.(newvalue_ismissingvalues, 1)
            ismissingvalues = getindex.(newvalue_ismissingvalues, 2)
            output[param] = AssocParam(param, newvalue, ismissingvalues, components, allcomponentsites, modelname, collect(paramsources[param]))
        elseif value isa Array{<:Any, 2}
            param ∈ asymmetric_pairparams && mirrormatrix!(value)
            newvalue, ismissingvalues = defaultmissing(value)
            if (!ignore_missingsingleparams 
                && !all([ismissingvalues[x,x] for x ∈ 1:size(ismissingvalues,1)])
                && any([ismissingvalues[x,x] for x ∈ 1:size(ismissingvalues,1)]))
                error("Partial missing values exist in diagonal of pair parameter ", param, ": ", [value[x,x] for x ∈ 1:size(ismissingvalues,1)], ".")
            end
            output[param] = PairParam(param, newvalue, ismissingvalues, components, allcomponentsites, modelname, collect(paramsources[param]))
        else
            error("Format for ", param, " is incorrect.")
        end
    end
    return output
end

function createparamarrays(filepaths::Array{String,1}, components::Array{String,1}, allcomponentsites::Array{Array{String,1},1}; verbose=false)
    # Returns Dict with all parameters in their respective arrays.
    checkfor_clashingheaders(filepaths)
    allparams = Dict{String,Any}()
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
        foundparams, paramtypes, sources = findparamsincsv(filepath, components, headerparams; verbose=verbose)
        foundcomponents = collect(keys(foundparams))
        foundparams = swapdictorder(foundparams)
        for headerparam ∈ headerparams
            !haskey(allparams, headerparam) && (allparams[headerparam] = createemptyparamsarray(paramtypes[headerparam], csvtype, components, allcomponentsites))
            !haskey(paramsources, headerparam) && (paramsources[headerparam] = Set{String}())
            if csvtype == singledata
                isempty(foundparams) && continue
                for (component, value) ∈ foundparams[headerparam]
                    currenttype = nonmissingtype(eltype(allparams[headerparam]))
                    if !(currenttype <: paramtypes[headerparam])
                        allparams[headerparam] = convert(Array{Union{Missing,paramtypes[headerparam]}}, allparams[headerparam])
                    end
                    idx = findfirst(isequal(component), components)
                    if allparams[headerparam] isa Array{<:Any,2}
                        !ismissing(allparams[headerparam][idx,idx]) && continue
                        allparams[headerparam][idx,idx] = value
                    else
                        !ismissing(allparams[headerparam][idx]) && continue
                        allparams[headerparam][idx] = value
                    end
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
                    if !(currenttype <: paramtypes[headerparam])
                        allparams[headerparam] = convert(Array{Union{Missing,paramtypes[headerparam]}}, allparams[headerparam])
                    end
                    idx1 = findfirst(isequal(componentpair[1]), components)
                    idx2 = findfirst(isequal(componentpair[2]), components)
                    !ismissing(allparams[headerparam][idx1,idx2]) && continue
                    allparams[headerparam][idx1,idx2] = value
                    !ismissing(sources[componentpair]) && push!(paramsources[headerparam], sources[componentpair])
                end
            end
            if csvtype == assocdata
                isempty(foundparams) && continue
                for (assocpair, value) ∈ foundparams[headerparam]
                    currenttype = nonmissingtype(eltype(first(allparams[headerparam])))
                    if !(currenttype <: paramtypes[headerparam])
                        allparams[headerparam] = convert(Array{Array{Union{Missing,paramtypes[headerparam]}}}, allparams[headerparam])
                    end
                    idx1 = findfirst(isequal(assocpair[1][1]), components)
                    idx2 = findfirst(isequal(assocpair[1][2]), components)
                    idx21 = findfirst(isequal(assocpair[2][1]), allcomponentsites[idx1])
                    idx22 = findfirst(isequal(assocpair[2][2]), allcomponentsites[idx2])
                    !ismissing(allparams[headerparam][idx1,idx2][idx21,idx22]) && continue
                    allparams[headerparam][idx1,idx2][idx21,idx22] = value
                    !ismissing(sources[assocpair]) && push!(paramsources[headerparam], sources[assocpair])
                end
            end
        end
    end
    return allparams, paramsources
end

function defaultmissing(array::Array)
    # Returns a non-missing Array with missing values replaced by default values.
    # Also returns an Array of Bools to retain overridden information.
    arraycopy = deepcopy(array)
    type = nonmissingtype(eltype(array))
    ismissingvalues = ismissing.(array)
    if type <: AbstractString
        arraycopy[ismissingvalues] .= ""
    elseif type <: Number
        arraycopy[ismissingvalues] .= 0
    else
        error("Unsupported type.")
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

function findparamsincsv(filepath::String, components::Array{String,1}, headerparams::Array{String,1}; columnreference="species", sitecolumnreference="site", sourcecolumnreference="source", verbose=false, ignore_missingsingleparams=false)
    # Returns a Dict with all matches in a particular file for one parameter.
    normalised_columnreference = normalisestring(columnreference)
    csvtype = readcsvtype(filepath)
    df = CSV.File(filepath; header=3)
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
    return lowercase(replace(str, ' ' => ""))
end

function findgroupsincsv(filepath::String, components::Array{String,1}; columnreference="species", groupcolumnreference="groups", verbose=false)
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
            
function readcsvheader(filepath::String; headerline = 3)
    # Returns array of filtered header strings at line 3.
    headers = split(getline(filepath, headerline), ',')
    ignorelist = ["source", "species", "dipprnumber", "smiles", "site"]
    return String.(filter(x -> replace.(lowercase(x), r"[ \d]" => "") ∉ ignorelist, headers))
end

function checkfor_clashingheaders(filepaths::Array{String,1})
    # Raises an error if the header of any assoc parameter clashes with a non-assoc parameter
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

function findsitesincsvs(filepaths::Array{String,1}, components::Array{String,1}; columnreference="species", sitecolumnreference="site", verbose=false)
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

function convertsingletopair(params::Array{T,1}; outputmissing=false) where T
    # Returns a missing square matrix with its diagonal matrix replaced by the given parameters. 
    paramslength = length(params)
    if outputmissing
        output = Array{Union{Missing,T}}(undef, paramslength, paramslength) .= missing
    else
        if T <: AbstractString
            output = Array{T}(undef, paramslength, paramslength) .= ""
        elseif T <: Number
            output = Array{T}(undef, paramslength, paramslength) .= 0
        else
            error("Data type ", T, " not supported.")
        end
    end
    for i = 1:paramslength
        output[i,i] = params[i]
    end
    return output
end

function mirrormatrix!(matrix::Array{<:Any,2})
    # Mirrors a square matrix.
    matrixsize = size(matrix)
    matrixsize[1] != matrixsize[2] && error("Matrix is not square.")
    for i ∈ 2:matrixsize[1], j ∈ 1:i-1
        lowervalue = matrix[i,j]
        uppervalue = matrix[j,i]
        if !ismissing(lowervalue) && !ismissing(uppervalue) && lowervalue != uppervalue
            error("Dissimilar non-zero entries exist across diagonal.")
        end
        !ismissing(lowervalue) && (matrix[j,i] = lowervalue)
        !ismissing(uppervalue) && (matrix[i,j] = uppervalue)
    end
    return matrix
end


function buildspecies(gccomponents::Array{<:Any,1}, models::Array{String,1}=String[]; usermodels::Array{String,1}=String[], modelname="", verbose=false)
    # The format for gccomponents is an arary of either the species name (if it
    # available in the OpenSAFT database, or a tuple consisting of the species
    # name, followed by a list of group => multiplicity pairs.  For example:
    # gccomponents = ["ethane",
    #                ("hexane", ["CH3" => 2, "CH2" => 4]),
    #                ("octane", ["CH3" => 2, "CH2" => 6])]
    BuildSpeciesType = Union{Tuple{String, Array{Pair{String, Int64},1}}, String, Tuple{String}}
    any(.!(isa.(gccomponents, BuildSpeciesType))) && error("The format of the components is incorrect.")
    filepaths = string.(vcat([(getdatabasepaths.(models)...)...], [(getuserpaths.(usermodels)...)...]))
    components = String[]
    allcomponentgroups = Array{Array{String,1},1}(undef, 0)
    allncomponentgroups = Array{Array{Int,1},1}(undef, 0)
    componentstolookup = String[]
    append!(componentstolookup, [x for x ∈ gccomponents[isa.(gccomponents, String)]])
    append!(componentstolookup, [first(x) for x ∈ gccomponents[isa.(gccomponents, Tuple{String})]])
    allfoundcomponentgroups = Dict{String,String}()
    for filepath in filepaths
        csvtype = readcsvtype(filepath)
        if csvtype != groupdata
            verbose && println("Skipping ", csvtype, " csv at ", filepath)
            continue
        end
        verbose && println("Searching for groups for components ", componentstolookup, " at ", filepath, "...")
        merge!(allfoundcomponentgroups, findgroupsincsv(filepath, componentstolookup; verbose=verbose))
    end
    for gccomponent ∈ gccomponents
        if gccomponent isa Tuple{String, Array{Pair{String, Int64},1}}
            append!(components, [gccomponent[1]])
            groupsandngroups = gccomponent[2]
        elseif gccomponent isa String
            !haskey(allfoundcomponentgroups, gccomponent) && error(gccomponent, " not found in any input csvs.")
            append!(components, [gccomponent])
            groupsandngroups = eval(Meta.parse(allfoundcomponentgroups[gccomponent]))
        elseif gccomponent isa Union{String, Tuple{String}}
            !haskey(allfoundcomponentgroups, first(gccomponent)) && error(gccomponent, " not found in any input csvs.")
            append!(components, [first(gccomponent)])
            groupsandngroups = eval(Meta.parse(allfoundcomponentgroups[first(gccomponent)]))
        end
        componentgroups = String[]
        ncomponentgroups = Int[]
        for (componentgroup, ncomponentgroup) ∈ groupsandngroups
            append!(componentgroups, [componentgroup])
            append!(ncomponentgroups, [ncomponentgroup])
        end
        append!(allcomponentgroups, [componentgroups])
        append!(allncomponentgroups, [ncomponentgroups])
    end
    if modelname == ""
        modelname = getmodelname(models, usermodels)
    end
    flattenedgroups = unique([(allcomponentgroups...)...])
    return GCParam(components, allcomponentgroups, allncomponentgroups, flattenedgroups, modelname)
end
