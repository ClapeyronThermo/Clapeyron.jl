@enum CSVType singledata pairdata assocdata groupdata

"""
    getfileextension(filepath)

A quick helper to get the file extension of any given path (without the dot).

# Examples
```julia-repl
julia> getfileextension("~/Desktop/text.txt")
"txt"
```
""" 
function getfileextension(filepath::AbstractString)
    path,ext = splitext(filepath)
    return String(chop(ext,head=1,tail=0))
end

"""
    getpaths(location; relativetodatabase=false)

Returns database paths that is optionally relative to Clapeyron.jl directory.
If path is a file, then return an Array containing a single path to that file.
If path is a directory, then return an Array containing paths to all csv files in that directory.

# Examples
```julia-repl
julia> getpaths("SAFT/PCSAFT"; relativetodatabase=true)
3-element Array{String,1}:
 "/home/user/.julia/packages/Clapeyron.jl/xxxxx/database/SAFT/PCSAFT/data_PCSAFT_assoc.csv"
 "/home/user/.julia/packages/Clapeyron.jl/xxxxx/database/SAFT/PCSAFT/data_PCSAFT_like.csv"
 "/home/user/.julia/packages/Clapeyron.jl/xxxxx/database/SAFT/PCSAFT/data_PCSAFT_unlike.csv"

```
"""
function getpaths(location::AbstractString; relativetodatabase::Bool=false)
    # We do not use realpath here directly because we want to make the .csv suffix optional.
    filepath = relativetodatabase ? normpath(dirname(pathof(Clapeyron)), "..", "database", location) : location
    isfile(filepath) && return [realpath(filepath)]
    isfile(filepath * ".csv") && return [realpath(filepath * ".csv")]
    if !isdir(filepath)
        relativetodatabase ? error("The path ", location, " does not exist in the Clapeyron database.") :
            error("The path ", location, " does not exist.")
    end
    files = joinpath.(filepath, readdir(filepath))
    return realpath.(files[isfile.(files) .& (getfileextension.(files) .== "csv")])
end

function flattenfilepaths(locations,userlocations)
    res = 
                vcat(
                    reduce(vcat,getpaths.(locations; relativetodatabase=true),init = String[]),
                    reduce(vcat,getpaths.(userlocations),init = String[]),
                    String[]
                    )
    return res
end

function getparams(components::Array{String,1}, locations::Array{String,1}=String[]; userlocations::Array{String,1}=String[], asymmetricparams::Array{String,1}=String[], ignore_missingsingleparams::Bool=false, verbose::Bool=false)
    # Gets all parameters from database.
    # locations is a list of paths relative to the Clapeyron database directory.
    # userlocations is a list of paths input by the user.
    # If parameters exist in multiple files, Clapeyron gives priority to files in later paths.
    # asymmetricparams is a list of parameters for which matrix reflection is disabled.
    # ignore_missingsingleparams gives users the option to disable component existence check in single params.
    filepaths = flattenfilepaths(locations,userlocations)
    allcomponentsites = findsitesincsvs(components, filepaths; verbose=verbose)
    allparams, paramsourcecsvs, paramsources = createparamarrays(components, filepaths, allcomponentsites; verbose=verbose)
    return packageparams(allparams, components, allcomponentsites, paramsourcecsvs, paramsources; asymmetricparams=asymmetricparams, ignore_missingsingleparams=ignore_missingsingleparams)
end

function getparams(groups::GroupParam, locations::Array{String,1}=String[]; userlocations::Array{String,1}=String[], asymmetricparams::Array{String,1}=String[], ignore_missingsingleparams::Bool=false, verbose::Bool=false)
    return getparams(groups.flattenedgroups, locations; userlocations=userlocations, asymmetricparams=asymmetricparams, ignore_missingsingleparams=ignore_missingsingleparams, verbose=verbose)
end

function getparams(components::String, locations::Array{String,1}=String[]; userlocations::Array{String,1}=String[], asymmetricparams::Array{String,1}=String[], ignore_missingsingleparams::Bool=false, verbose::Bool=false)
    return getparams([components], locations; userlocations=userlocations, asymmetricparams=asymmetricparams, ignore_missingsingleparams=ignore_missingsingleparams, verbose=verbose)
end

function packageparams(allparams::Dict, components::Array{String,1}, allcomponentsites::Array{Array{String,1},1}, paramsourcecsvs::Dict{String,Set{String}}, paramsources::Dict{String,Set{String}}; asymmetricparams::Array{String,1}=String[], ignore_missingsingleparams::Bool=false)
    # Package params into their respective Structs.
    output = Dict{String,ClapeyronParam}()
    for (param, value) ∈ allparams
        if value isa Vector 
            newvalue, ismissingvalues = defaultmissing(value)
            if !ignore_missingsingleparams && any(ismissingvalues)
                error("Missing values exist in single parameter ", param, ": ", value, ".")
            end
            output[param] = SingleParam(param, newvalue, ismissingvalues, components, allcomponentsites, collect(paramsourcecsvs[param]), collect(paramsources[param]))
        elseif value isa Matrix{<:Array}
            param ∉ asymmetricparams && mirrormatrix!(value) 
            newvalue_ismissingvalues = defaultmissing.(value)
            newvalue = getindex.(newvalue_ismissingvalues, 1)
            ismissingvalues = getindex.(newvalue_ismissingvalues, 2)
            output[param] = AssocParam(param, newvalue, ismissingvalues, components, allcomponentsites, collect(paramsourcecsvs[param]), collect(paramsources[param]))
        elseif value isa Matrix
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

function param_type(t1,t2)
    t_promoted = promote_type(t1,t2)
    local res
    if t_promoted !== Any
        res =  Union{Missing,t_promoted}
    else
        res = Union{Missing,t1,t2}
    end
    return res
end

#try not to promote to Int64. this will error if a user pass integer values different from 0 and 1
function param_type(::Type{Bool},::Type{Int})
    return Union{Bool,Missing}
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
        headerparams = readheaderparams(filepath)
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
                    tx = param_type(paramtypes[headerparam],currenttype)
                    allparams[headerparam] = convert(Array{tx}, allparams[headerparam])
                    idx = findfirst(isequal(component), components)
                    if allparams[headerparam] isa Matrix
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
                if allparams[headerparam] isa Vector
                    allparams[headerparam] = convertsingletopair(allparams[headerparam], true)
                end
                isempty(foundparams) && continue
                for (componentpair, value) ∈ foundparams[headerparam]
                    currenttype = nonmissingtype(eltype(allparams[headerparam]))
                    tx = param_type(paramtypes[headerparam],currenttype)
                    allparams[headerparam] = convert(Matrix{tx}, allparams[headerparam])
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
                    tx = param_type(paramtypes[headerparam],currenttype)
                    allparams[headerparam] = convert(Array{Array{tx,2},2}, allparams[headerparam])
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

function defaultmissing(array::Array{<:Number},defaultvalue)
    return deepcopy(array),Array(ismissing.(array))
end

function defaultmissing(array::Array{<:AbstractString},defaultvalue)
    return string.(array),Array(ismissing.(array))
end

function defaultmissing(array::Array{String},defaultvalue)
    return deepcopy(array),Array(ismissing.(array))
end

function defaultmissing(array::Array{Union{Missing, T}},defaultvalue="") where T <:AbstractString
    return string.(coalesce.(array,defaultvalue)),Array(ismissing.(array))
end

function defaultmissing(array::Array{Union{Missing, T}},defaultvalue=zero(eltype(array))) where T<:Number
    return coalesce.(array,defaultvalue),Array(ismissing.(array))
end

function defaultmissing(array::Array{Union{Missing, Bool}},defaultvalue=zero(eltype(array)))
    return [coalesce(i,defaultvalue) for i in array],Array(ismissing.(array))
end

#if an array with only missings is passed, the Resulting ClapeyronParam will be 
#of the type that this function returns
function defaultmissing(array::Array{Missing},defaultvalue=0.0)
    return array,Array(ismissing.(array))
end

function defaultmissing(array::Array{Union{Missing,T1,T2}},defaultvalue="") where {T1 <:Number,T2<:AbstractString}
    return string.(coalesce.(array,defaultvalue)),Array(ismissing.(array))
end

function defaultmissing(array::Array{Any},defaultvalue="")
    return string.(coalesce.(array,defaultvalue)),Array(ismissing.(array))
end

function defaultmissing(array,defaultvalue="")
    throw("Unsupported array element type  $(typeof(array))")
end


function swapdictorder(dict)
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
function findparamsincsv(components::Array{String,1},
    filepath::AbstractString,
    headerparams::Array{String,1};
    columnreference::AbstractString="species",
    sitecolumnreference::AbstractString="site",
    sourcecolumnreference::AbstractString="source",
    verbose::Bool=false,
    normalisecomponents::Bool=true)
    # Returns a Dict with all matches in a particular file for one parameter.
    normalised_columnreference = normalisestring(columnreference)
    normalised_columnreference1 = normalised_columnreference * '1'
    normalised_columnreference2 = normalised_columnreference * '2'

    csvtype = readcsvtype(filepath)
    df = CSV.File(filepath; header=3, pool=0,lazystrings=true,silencewarnings=true)
    csvheaders = String.(Tables.columnnames(df))
    normalised_components = normalisestring.(components; isactivated=normalisecomponents)
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
            foundcomponentidx = findfirst(isequal(normalisestring(component; isactivated=normalisecomponents)), normalised_components)
            isnothing(foundcomponentidx) && continue
            verbose && print("Found component: ", component)
            component = components[foundcomponentidx]
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
        lookupcolumnindex1 = findfirst(isequal(normalised_columnreference1), normalised_csvheaders)
        lookupcolumnindex2 = findfirst(isequal(normalised_columnreference2), normalised_csvheaders)
        isnothing(lookupcolumnindex1) && error("Header ",normalised_columnreference1, " not found.")
        isnothing(lookupcolumnindex2) && error("Header ", normalised_columnreference2, " not found.")
        lookupcolumn1 = Symbol(csvheaders[lookupcolumnindex1])
        lookupcolumn2 = Symbol(csvheaders[lookupcolumnindex2])
        for row ∈ Tables.rows(df)
            component1 = row[lookupcolumn1]
            component2 = row[lookupcolumn2]
            foundcomponentidx1 = findfirst(isequal(normalisestring(component1; isactivated=normalisecomponents)), normalised_components)
            foundcomponentidx2 = findfirst(isequal(normalisestring(component2; isactivated=normalisecomponents)), normalised_components)
            (isnothing(foundcomponentidx1) || isnothing(foundcomponentidx2)) && continue
            verbose && print("Found component pair: ", (component1, component2))
            componentpair = (components[foundcomponentidx1], components[foundcomponentidx2])
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
        lookupcolumnindex1 = findfirst(isequal(normalised_columnreference1), normalised_csvheaders)
        lookupcolumnindex2 = findfirst(isequal(normalised_columnreference2), normalised_csvheaders)
        isnothing(lookupcolumnindex1) && error("Header ", normalised_columnreference1, " not found.")
        isnothing(lookupcolumnindex2) && error("Header ", normalised_columnreference2, " not found.")
        lookupcolumn1 = Symbol(csvheaders[lookupcolumnindex1])
        lookupcolumn2 = Symbol(csvheaders[lookupcolumnindex2])
        normalised_sitecolumnreference = normalisestring(sitecolumnreference)
        normalised_sitecolumnreference1 = normalised_sitecolumnreference * '1'
        normalised_sitecolumnreference2 = normalised_sitecolumnreference * '2'
        lookupsitecolumnindex1 = findfirst(isequal(normalised_sitecolumnreference1), normalised_csvheaders)
        lookupsitecolumnindex2 = findfirst(isequal(normalised_sitecolumnreference2), normalised_csvheaders)
        isnothing(lookupcolumnindex1) && error("Header ", normalised_sitecolumnreference1, " not found.")
        isnothing(lookupcolumnindex2) && error("Header ", normalised_sitecolumnreference2, " not found.")
        lookupsitecolumn1 = Symbol(csvheaders[lookupsitecolumnindex1])
        lookupsitecolumn2 = Symbol(csvheaders[lookupsitecolumnindex2])
        for row ∈ Tables.rows(df)
            component1 = row[lookupcolumn1]
            component2 = row[lookupcolumn2]
            foundcomponentidx1 = findfirst(isequal(normalisestring(component1; isactivated=normalisecomponents)), normalised_components)
            foundcomponentidx2 = findfirst(isequal(normalisestring(component2; isactivated=normalisecomponents)), normalised_components)
            (isnothing(foundcomponentidx1) || isnothing(foundcomponentidx2)) && continue
            site1 = row[lookupsitecolumn1]
            site2 = row[lookupsitecolumn2]
            verbose && print("Found assoc pair: ", ((component1, component2),(site1, site2)))
            assocpair = ((components[foundcomponentidx1], components[foundcomponentidx2]), (site1, site2))
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

function normalisestring(str::AbstractString; isactivated::Bool=true, tofilter::Regex=r"[ \-\_]", changecase::Bool=true)
    !isactivated && return str
    normalisedstring = replace(str, tofilter => "")
    changecase && (normalisedstring = lowercase(normalisedstring))
    return normalisedstring
end

function readcsvtype(filepath::AbstractString)
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

function getline(filepath::AbstractString, selectedline::Int)
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

const readheaderparams_ignorelist = ["source", "species", "dipprnumber", "smiles", "site"]

function readheaderparams(filepath::AbstractString; headerline::Int = 3)
    # Returns array of filtered header strings at line 3.
    ignorelist = readheaderparams_ignorelist
    headers = split(getline(filepath, headerline), ',')
    if last(headers) == ""
        pop!(headers)
    end
    return String.(filter(x -> normalisestring(x; tofilter=r"[ \-\_\d]") ∉ ignorelist, headers))
end

function checkfor_clashingheaders(filepaths::Array{String,1})
    # Raises an error if the header of any assoc parameter clashes with a non-assoc parameter.
    headerparams = String[]
    headerparams_assoc = String[]
    for filepath in filepaths
        csvtype = readcsvtype(filepath)
        if csvtype == singledata || csvtype == pairdata
            append!(headerparams, readheaderparams(filepath))
        elseif csvtype == assocdata
            append!(headerparams_assoc, readheaderparams(filepath))
        end
    end
    clashingheaders = intersect(headerparams, headerparams_assoc)
    !isempty(clashingheaders) && error("Headers ", clashingheaders, " appear in both loaded assoc and non-assoc files.")
end

function findsitesincsvs(components::Array{String,1}, filepaths::Array{String,1}; columnreference::AbstractString="species", sitecolumnreference::AbstractString="site", verbose::Bool=false, normalisecomponents::Bool=true)
    # Look for all relevant sites in the database.
    # Note that this might not necessarily include all sites associated with a component.
    normalised_components = normalisestring.(components; isactivated=normalisecomponents)
    normalised_columnreference = normalisestring(columnreference)
    normalised_columnreference1 = normalised_columnreference * '1'
    normalised_columnreference2 = normalised_columnreference * '2'
    normalised_sitecolumnreference = normalisestring(sitecolumnreference)
    sites = Dict(components .=> [Set{String}() for _ ∈ 1:length(components)])
    for filepath ∈ filepaths
        csvtype = readcsvtype(filepath)
        csvtype != assocdata && continue
        df = CSV.File(filepath; header=3,lazystrings=true)
        csvheaders = String.(Tables.columnnames(df))
        normalised_csvheaders = normalisestring.(String.(Tables.columnnames(df)))
        lookupcolumnindex1 = findfirst(isequal(normalised_columnreference1), normalised_csvheaders)
        lookupcolumnindex2 = findfirst(isequal(normalised_columnreference2), normalised_csvheaders)
        isnothing(lookupcolumnindex1) && error("Header ", normalised_columnreference1, " not found.")
        isnothing(lookupcolumnindex2) && error("Header ", normalised_columnreference2, " not found.")
        lookupcolumn1 = Symbol(csvheaders[lookupcolumnindex1])
        lookupcolumn2 = Symbol(csvheaders[lookupcolumnindex2])

        normalised_sitecolumnreference = normalisestring(sitecolumnreference)
        normalised_sitecolumnreference1 = normalised_sitecolumnreference * '1'
        normalised_sitecolumnreference2 = normalised_sitecolumnreference * '2'

        lookupsitecolumnindex1 = findfirst(isequal(normalised_sitecolumnreference1), normalised_csvheaders)
        lookupsitecolumnindex2 = findfirst(isequal(normalised_sitecolumnreference2), normalised_csvheaders)
        isnothing(lookupcolumnindex1) && error("Header ", normalised_sitecolumnreference1, " not found.")
        isnothing(lookupcolumnindex2) && error("Header ", normalised_sitecolumnreference2, " not found.")
        lookupsitecolumn1 = Symbol(csvheaders[lookupsitecolumnindex1])
        lookupsitecolumn2 = Symbol(csvheaders[lookupsitecolumnindex2])

        for row ∈ Tables.rows(df)
            component1 = row[Symbol(lookupcolumn1)]
            component2 = row[Symbol(lookupcolumn2)]
            foundcomponentidx1 = findfirst(isequal(normalisestring(component1; isactivated=normalisecomponents)), normalised_components)
            foundcomponentidx2 = findfirst(isequal(normalisestring(component2; isactivated=normalisecomponents)), normalised_components)
            (isnothing(foundcomponentidx1) || isnothing(foundcomponentidx2)) && continue
            push!(sites[components[foundcomponentidx1]], row[Symbol(lookupsitecolumn1)])
            push!(sites[components[foundcomponentidx2]], row[Symbol(lookupsitecolumn2)])
        end
    end
    output = Array{Array{String,1}}(undef, 0)
    for component ∈ components
        push!(output, collect(sites[component]))
    end
    verbose && println("Found sites for ", components, " are ", output, ".")
    return output
end


function findgroupsincsv(components::Array{String,1},
    filepath::AbstractString; 
    columnreference::AbstractString="species", 
    groupcolumnreference::AbstractString="groups", 
    verbose::Bool=false,
    normalisecomponents=true)

    # Returns a Dict with the group string that will be parsed in buildspecies.
    csvtype = readcsvtype(filepath)
    csvtype != groupdata && return Dict{String,String}()
    normalised_components = normalisestring.(components; isactivated=normalisecomponents)
    normalised_columnreference = normalisestring(columnreference)
    normalised_groupcolumnreference = normalisestring(groupcolumnreference)
    csvtype = readcsvtype(filepath)
    df = CSV.File(filepath; header=3,lazystrings=true)
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
    foundcomponentidx = findfirst(isequal(normalisestring(component; isactivated=normalisecomponents)), normalised_components)
    isnothing(foundcomponentidx) && continue
    verbose && print("Found component: ", component)
    component = components[foundcomponentidx]
    foundgroups[component] = row[lookupgroupcolumn]
    verbose && println(" with groups ", foundgroups[component])
    end
    return foundgroups
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
 
"""
    convertsingletopair(params::Vector,outputmissing=false)

For numbers, this is equal to Matrix(Diagonal(params)).
For strings, this is equal to a matrix filled with "", whose diagonal is `params`

"""
function convertsingletopair(params::Array{T,1}, outputmissing::Bool=false) where T
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

function mirrormatrix!(matrix::Matrix{T}) where T
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
function GroupParam(gccomponents, grouplocations::Array{String,1}=String[]; usergrouplocations::Array{String,1}=String[], verbose::Bool=false)
    if (!(gccomponents isa PARSED_GROUP_VECTOR_TYPE)) | iszero(length(grouplocations)) | iszero(length(usergrouplocations)) | verbose
        return _GroupParam(gccomponents,grouplocations;usergrouplocations,verbose)
    else
         return GroupParam(gccomponents)
    end        
end
    
function _GroupParam(gccomponents, grouplocations::Array{String,1}=String[]; usergrouplocations::Array{String,1}=String[], verbose::Bool=false)
    # The format for gccomponents is an arary of either the species name (if it
    # available in the Clapeyron database, or a tuple consisting of the species
    # name, followed by a list of group => multiplicity pairs.  For example:
    # gccomponents = ["ethane",
    #                ("hexane", ["CH3" => 2, "CH2" => 4]),
    #                ("octane", ["CH3" => 2, "CH2" => 6])]
    BuildSpeciesType = Union{Tuple{String, Array{Pair{String, Int64},1}}, String, Tuple{String}}
    any(.!(isa.(gccomponents, BuildSpeciesType))) && error("The format of the components is incorrect.")
    filepaths = flattenfilepaths(grouplocations,usergrouplocations)
    componentstolookup = String[]
    #componentstolookup = filter(x-> x isa Union{String,Tuple{String}},gccomponents)
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
    gccomponents_parsed = PARSED_GROUP_VECTOR_TYPE(undef,length(gccomponents))
    for (i,gccomponent) ∈ pairs(gccomponents)
        if gccomponent isa Tuple{String, Array{Pair{String, Int64},1}}
            gccomponents_parsed[i] = gccomponent
        elseif gccomponent isa String
            !haskey(allfoundcomponentgroups, gccomponent) && error("Predefined component ", gccomponent, " not found in any group input csvs.")
            groupsandngroups = eval(Meta.parse(allfoundcomponentgroups[gccomponent]))
            gccomponents_parsed[i] = (gccomponent,groupsandngroups)
        elseif gccomponent isa Tuple{String}
            !haskey(allfoundcomponentgroups, first(gccomponent)) && error("Predefined component ", gccomponent, " not found in any group input csvs.")
            groupsandngroups = eval(Meta.parse(allfoundcomponentgroups[first(gccomponent)]))
            gccomponents_parsed[i] = (first(gccomponent),groupsandngroups)
        end
    end
    return GroupParam(gccomponents_parsed,groupsourcecsvs)
end

export getparams
