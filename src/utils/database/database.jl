using CSV, Tables
include("OpenSAFTParams.jl")
include("checks.jl")
include("visualisation.jl")

@enum ParameterType singledata pairdata assocdata groupdata

function getfileextension(filepath::String)
    # Quick helper function to get the file extension of any given path.
    dotpos = findlast(isequal('.'), filepath)
    isnothing(dotpos) && return ""
    return filepath[dotpos+1:end]
end

function getdatabasepaths(model::String)
    # Returns database paths relative to OpenSAFT.jl directory.
    # If path is a file, then return an Array containing a single path to that file.
    # If path is a directory, then return an Array containing paths to all csv files in that directory.
    path = joinpath(dirname(pathof(OpenSAFT)), "../database", model)
    isfile(path) && return [path]
    isfile(path * ".csv") && return [path * ".csv"]
    !isdir(path) && error("The directory ", model, " does not exist in the OpenSAFT database.")
    files = joinpath.(path, readdir(path))
    return files[isfile.(files) .& (getfileextension.(files) .== "csv")]
end

function getuserpaths(model::String)
    # If path is a file, then return an Array containing a single path to that file.
    # If path is a directory, then return an Array containing paths to all csv files in that directory.
    path = model
    isfile(path) && return [path]
    isfile(path * ".csv") && return [path * ".csv"]
    !isdir(path) && error("The directory ", path, " does not exist.")
    files = joinpath.(path, readdir(path))
    return files[isfile.(files) .& (getfileextension.(files) .== "csv")]
end

function getmodelname(models::Array{String,1}, usermodels::Array{String,1})
    # Try to guess the name of the model.
    # It will take the name of the first given directory, checking models before usermodels.
    if !isempty(models)
        for model in models
            path = joinpath(dirname(pathof(OpenSAFT)), "../database", model)
            ispath(path) && return basename(path)
        end
    end
    if !isempty(usermodels)
        for usermodel in usermodels
            path = usermodel
            ispath(path) && return basename(path)
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
    sites = findsites(filepaths, components)
    allparams, paramsources = findparams(filepaths, components, sites; verbose=verbose)
    if modelname == ""
        modelname = getmodelname(models, usermodels)
    end
    finaldict = packageparams(allparams, components, sites, paramsources, modelname; asymmetric_pairparams=asymmetric_pairparams, ignore_missingsingleparams=ignore_missingsingleparams)
    return finaldict
end

function packageparams(allparams::Dict, components::Array{String,1}, sites::Array{Array{String,1},1}, paramsources::Dict{String,Set{String}}, modelname::String; asymmetric_pairparams::Array{String,1}=String[], ignore_missingsingleparams=false)
    # Package params into their respective Structs.
    output = Dict()
    for (param, value) in allparams
        if typeof(value) <: Array{<:Any,1}
            newvalue, ismissingvalues = defaultmissing(value)
            if !ignore_missingsingleparams && any(ismissingvalues)
                error("Missing values exist in single parameter ", param, ": ", value, ".")
            end
            output[param] = SingleParams(param, newvalue, ismissingvalues, components, modelname, collect(paramsources[param]))
        elseif typeof(value) <: Array{<:Array,2}
            newvalue_ismissingvalues = defaultmissing.(value)
            newvalue = getindex.(newvalue_ismissingvalues, 1)
            ismissingvalues = getindex.(newvalue_ismissingvalues, 2)
            output[param] = AssocParams(param, newvalue, ismissingvalues, components, sites, modelname, collect(paramsources[param]))
        elseif typeof(value) <: Array{<:Any, 2}
            param in asymmetric_pairparams && mirrormatrix!(value)
            newvalue, ismissingvalues = defaultmissing(value)
            if (!ignore_missingsingleparams 
                && !all([ismissingvalues[x,x] for x in 1:size(ismissingvalues,1)])
                && any([ismissingvalues[x,x] for x in 1:size(ismissingvalues,1)]))
                error("Partial missing values exist in diagonal of pair parameter ", param, ": ", [value[x,x] for x in 1:size(ismissingvalues,1)], ".")
            end
            output[param] = PairParams(param, newvalue, ismissingvalues, components, modelname, collect(paramsources[param]))
        else
            error("Format for ", param, " is incorrect.")
        end
    end
    return output
end

function findparams(filepaths::Array{String,1}, components::Array{String,1}, sites::Array{Array{String,1},1}; verbose=false)
    # Returns Dict with all parameters in their respective arrays.
    checkfor_clashingheaders(filepaths)
    allparams = Dict{String,Any}()
    paramsources = Dict{String,Set{String}}()
    for filepath in filepaths
        type = readtype(filepath)
        headerparams = readheader(filepath)
        foundparams, paramtypes, sources = searchfor(filepath, components, headerparams; verbose=verbose)
        isempty(foundparams) && continue
        foundcomponents = collect(keys(foundparams))
        foundparams = swapdictorder(foundparams)
        for headerparam in headerparams
            if !haskey(allparams, headerparam)
                allparams[headerparam] = createemptyparamsarray(paramtypes[headerparam], type, components, sites)
            end
            if type == singledata
                for (component, value) in foundparams[headerparam]
                    currenttype = nonmissingtype(eltype(allparams[headerparam]))
                    if !(currenttype <: paramtypes[headerparam])
                        allparams[headerparam] = convert(Array{Union{Missing,paramtypes[headerparam]}}, allparams[headerparam])
                    end
                    idx = findfirst(isequal(component), components)
                    if typeof(allparams[headerparam]) <: Array{<:Any,2}
                        allparams[headerparam][idx,idx] = value
                    else
                        allparams[headerparam][idx] = value
                    end
                    !haskey(paramsources, headerparam) && (paramsources[headerparam] = Set())
                    !ismissing(sources[component]) && push!(paramsources[headerparam], sources[component])
                end
            end
            if type == pairdata
                if typeof(allparams[headerparam]) <: Array{<:Any,1}
                    allparams[headerparam] = convertsingletopair(allparams[headerparam])
                end
                for (componentpair, value) in foundparams[headerparam]
                    currenttype = nonmissingtype(eltype(allparams[headerparam]))
                    if !(currenttype <: paramtypes[headerparam])
                        allparams[headerparam] = convert(Array{Union{Missing,paramtypes[headerparam]}}, allparams[headerparam])
                    end
                    idx1 = findfirst(isequal(componentpair[1]), components)
                    idx2 = findfirst(isequal(componentpair[2]), components)
                    allparams[headerparam][idx1,idx2] = value
                    !haskey(paramsources, headerparam) && (paramsources[headerparam] = Set())
                    !ismissing(sources[componentpair]) && push!(paramsources[headerparam], sources[componentpair])
                end
            end
            if type == assocdata
                for (assocpair, value) in foundparams[headerparam]
                    currenttype = nonmissingtype(eltype(allparams[headerparam][1,1]))
                    if !(currenttype <: paramtypes[headerparam])
                        allparams[headerparam] = convert(Array{Array{Union{Missing,paramtypes[headerparam]}}}, allparams[headerparam])
                    end
                    idx1 = findfirst(isequal(assocpair[1][1]), components)
                    idx2 = findfirst(isequal(assocpair[1][2]), components)
                    idx21 = findfirst(isequal(assocpair[2][1]), sites[idx1])
                    idx22 = findfirst(isequal(assocpair[2][2]), sites[idx2])
                    allparams[headerparam][idx1,idx2][idx21,idx22] = value
                    !haskey(paramsources, headerparam) && (paramsources[headerparam] = Set())
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
    # Swap the first two level in a nested dictionary.
    # Note that there is no checking done to ensure that Dict format is correct
    isempty(dict) && return dict
    output = Dict()
    outerkeys = keys(dict)
    innerkeys = keys(dict[collect(outerkeys)[1]])
    for innerkey in innerkeys, outerkey in outerkeys
        if !haskey(output, innerkey)
            output[innerkey] = Dict{Any,Any}(outerkey => dict[outerkey][innerkey])
        end
        push!(output[innerkey], outerkey => dict[outerkey][innerkey])
    end
    return output
end

function searchfor(filepath::String, components::Array{String,1}, headerparams::Array{String,1}; columnreference="species", sitecolumnreference="site", sourcecolumnreference="source", verbose=false, ignore_missingsingleparams=false)
    # Returns a Dict with all matches in a particular file for one parameter.
    normalised_columnreference = lowercase(replace(columnreference, ' ' => ""))
    type = readtype(filepath)
    verbose && println("Searching for ", string(type), " headers ", headerparams, " for components ", components, " at ", filepath, "...")
    df = CSV.File(filepath; header=3)
    csvheaders = replace.(lowercase.(String.(Tables.columnnames(df))), ' ' => "")

    foundvalues = Dict()
    paramtypes = Dict(headerparams .=> [Tables.columntype(df, Symbol(x)) for x in headerparams])
    sources = Dict()

    normalised_sourcecolumnreference = lowercase(replace(sourcecolumnreference, ' ' => ""))
    getsources = false
    if normalised_sourcecolumnreference in csvheaders
        getsources = true
        sourcecolumn = Symbol(csvheaders[findfirst(isequal(normalised_sourcecolumnreference), csvheaders)])
    end

    if type == singledata
        lookupcolumn = Symbol(csvheaders[findfirst(isequal(normalised_columnreference), csvheaders)])
        for row in Tables.rows(df)
            component = row[lookupcolumn]
            if component in components
                verbose && print("Found component: ", component)
                foundvalues[component] = Dict()
                for headerparam in headerparams
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
        end
    elseif type == pairdata
        lookupcolumn1 = Symbol(csvheaders[findfirst(x -> x[1:end-1] == normalised_columnreference && x[end] == '1', csvheaders)])
        lookupcolumn2 = Symbol(csvheaders[findfirst(x -> x[1:end-1] == normalised_columnreference && x[end] == '2', csvheaders)])
        for row in Tables.rows(df)
            component1 = row[lookupcolumn1]
            component2 = row[lookupcolumn2]
            if component1 in components && component2 in components
                componentpair = (component1, component2)
                verbose && print("Found component pair: ", componentpair)
                foundvalues[componentpair] = Dict()
                for headerparam in headerparams
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
        end
    elseif type == assocdata
        normalised_sitecolumnreference = lowercase(replace(sitecolumnreference, ' ' => ""))
        lookupcolumn1 = Symbol(csvheaders[findfirst(x -> x[1:end-1] == normalised_columnreference && x[end] == '1', csvheaders)])
        lookupcolumn2 = Symbol(csvheaders[findfirst(x -> x[1:end-1] == normalised_columnreference && x[end] == '2', csvheaders)])
        lookupsitecolumn1 = Symbol(csvheaders[findfirst(x -> x[1:end-1] == normalised_sitecolumnreference && x[end] == '1', csvheaders)])
        lookupsitecolumn2 = Symbol(csvheaders[findfirst(x -> x[1:end-1] == normalised_sitecolumnreference && x[end] == '2', csvheaders)])
        for row in Tables.rows(df)
            component1 = row[lookupcolumn1]
            component2 = row[lookupcolumn2]
            if component1 in components && component2 in components
                site1 = row[lookupsitecolumn1]
                site2 = row[lookupsitecolumn2]
                assocpair = ((component1, component2), (site1, site2))
                verbose && print("Found assoc pair: ", assocpair)
                foundvalues[assocpair] = Dict()
                for headerparam in headerparams
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
        end
    else
        error("File is of type ", String(type), " and cannot be read with this function.")
    end
    return foundvalues, paramtypes, sources
end

function readtype(filepath::String)
    # Searches for type from second line of CSV.
    keywords = ["like", "single", "unlike", "pair", "assoc", "group"]
    words = split(lowercase(rstrip(getline(filepath, 2), ',')), ' ')
    foundkeywords = intersect(words, keywords)
    length(foundkeywords) > 1 && error("Multiple keywords found in database ", filepath, ": ", foundkeywords)
    isempty(foundkeywords) && error("Unable to determine type of database", filepath, ". Check that keyword is present on Line 2.")
    foundkeywords[1] == "single" && return singledata
    foundkeywords[1] == "like" && return singledata
    foundkeywords[1] == "pair" && return pairdata
    foundkeywords[1] == "unlike" && return pairdata
    foundkeywords[1] == "assoc" && return assocdata
    foundkeywords[1] == "group" && return groupdata
end

function getline(filepath::String, selectedline::Int)
    # Simple function to return text from filepath at selectedline.
    open(filepath) do file
        linecount = 1
        for line in eachline(file)
            linecount == selectedline && return line
            linecount += 1
        end
        error("Selected line number exceeds number of lines in file")
    end
end
            
function readheader(filepath::String; headerline = 3)
    # Returns array of filtered header strings at line 3.
    headers = split(getline(filepath, headerline), ',')
    ignorelist = ["source", "species", "dipprnumber", "smiles", "site"]
    return String.(filter(x -> replace.(lowercase(x), r"[ \d]" => "") ∉ ignorelist, headers))
end

function findsites(filepaths::Array{String,1}, components::Array{String,1}; columnreference="species", sitecolumnreference="site", verbose=false)
    normalised_columnreference = lowercase(replace(columnreference, ' ' => ""))
    normalised_sitecolumnreference = lowercase(replace(sitecolumnreference, ' ' => ""))
    sites = Dict(components .=> [Set() for x in 1:length(components)])
    for filepath in filepaths
        type = readtype(filepath)
        type != assocdata && continue
        df = CSV.File(filepath; header=3)
        csvheaders = replace.(lowercase.(String.(Tables.columnnames(df))), ' ' => "")
        lookupcolumn1 = csvheaders[findfirst(x -> x[1:end-1] == normalised_columnreference && x[end] == '1', csvheaders)]
        lookupcolumn2 = csvheaders[findfirst(x -> x[1:end-1] == normalised_columnreference && x[end] == '2', csvheaders)]
        lookupsitecolumn1 = csvheaders[findfirst(x -> x[1:end-1] == normalised_sitecolumnreference && x[end] == '1', csvheaders)]
        lookupsitecolumn2 = csvheaders[findfirst(x -> x[1:end-1] == normalised_sitecolumnreference && x[end] == '2', csvheaders)]
        for row in Tables.rows(df)
            component1 = row[Symbol(lookupcolumn1)]
            component2 = row[Symbol(lookupcolumn2)]
            if component1 in components && component2 in components
                push!(sites[component1], row[Symbol(lookupsitecolumn1)])
                push!(sites[component2], row[Symbol(lookupsitecolumn2)])
            end
        end
    end
    output = Array{Array{String,1}}(undef, 0)
    for component in components
        push!(output, collect(sites[component]))
    end
    verbose && println("Found sites for ", components, " are ", output, ".")
    return output
end

function createemptyparamsarray(datatype::Type, type::ParameterType, components::Array{String,1}, sites::Array{Array{String,1},1})
    # Creates a missing array of the appropriate size.
    componentslength = length(components)
    type == singledata && return (Array{Union{Missing,datatype}}(undef, componentslength) .= missing)
    type == pairdata && return (Array{Union{Missing,datatype}}(undef, componentslength, componentslength) .= missing)
    if type == assocdata
        output = Array{Array{Union{Missing,datatype},2},2}(undef, componentslength, componentslength)
        for i in 1:componentslength, j in 1:componentslength
            output[i,j] = (Array{Union{Missing,datatype}}(undef, length(sites[i]), length(sites[j])) .= missing)
        end
        return output
    end
end

function createemptyparamsarray(type::ParameterType, components::Array{String,1}, sites::Array{Array{String,1},1})
    return createemptyparamsarray(Any, type, components, sites)
end

function convertsingletopair(params::Array{T,1}) where T
    # Returns a diagonal matrix with the given parameters. 
    paramslength = length(params)
    output = Array{Union{Missing,T}}(undef, paramslength, paramslength) .= missing
    for i = 1:paramslength
        output[i,i] = params[i]
    end
    return output
end

function mirrormatrix!(matrix::Array{<:Any,2})
    # Mirrors a square matrix.
    matrixsize = size(matrix)
    matrixsize[1] != matrixsize[2] && error("Matrix is not square.")
    for i in 2:matrixsize[1], j in 1:i-1
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
