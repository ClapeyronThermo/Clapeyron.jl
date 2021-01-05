@enum ParameterType likedata pairdata assocdata groupdata

abstract type OpenSAFTParams end
struct LikeParams(T) <: OpenSAFTParams
    name::String
    model::String
    values::Array{T,1}
    components::Array{String,1}
end

struct PairParams(T) <: OpenSAFTParams
    name::String
    model::String
    values::Array{T,2}
    components::Array{String,1}
end

struct AssocParams(T) <: OpenSAFTParams
    name::String
    model::String
    values::Array{T,1}
    components::Array{String,1}
    sites::Array{Array{String,1},1}
end

struct GroupParams <: OpenSAFTParams
    model::String
    components::Array{String,1}
    groups::Array{Array{String,1},1}
    groupmultiplicities::Array{Array{Int64,1},1}
end

function checkfor_clashingheaders(locations::Array{String,1})
    # Raises an error if the header of any assoc parameter clashes with a non-assoc parameter
    headerparams_like = Array{Array{String,1},1}
    headerparams_pair = Array{Array{String,1},1}
    headerparams_assoc = Array{Array{String,1},1}
    for location in locations
        type = readtype(location)
        if type == likedata || type == pairdata
            append!(headerparams, readheader(location))
        if type == assocdata
            append!(hearderparams_assoc, readheader(location))
        end
    end
    clashingheaders = intersect(headerparams_like, headerparams_pair, headerparams_assoc)
    isempty(clashingheaders) && error("Headers ", clashingheaders, " appear in both loaded asssoc and non-assoc files.")
end

function checkfor_likecompleteness(locations::Array{String,1}, components::Array{String,1})
    # Raises an error if any component is not present in any like databases.
end

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
    !isdir(directory) && error("The directory ", model, " does not exist in the OpenSAFT database.")
    files = joinpath.(directory, readdir(directory))
    return files[isfile.(files) .& (getfileextension.(files) .== "csv")]
end

function getparams(components::Array{String,1}, models::Array{String,1}=[]; filepaths::Array{String,1}=[], asymmetric_pairparams::Array{String,1}=[], ignore_missinglikeparams=false)
    # Gets all parameters from database.
    # models is a list of paths relative to the OpenSAFT database directory.
    # filepaths is a list of paths input by the user.
    # If parameters exist in multiple files, OpenSAFT gives priority to files in later paths.
    # asymmetric_pairparms is a list of parameters with array reflection turned off.
    # error_onmissingparams gives users the option to disable component existence check in like params.
    locations = vcat([(getdatabasepaths.(models)...)...], filepaths)
    ignore_missinglikeparams && checkfor_likecompleteness(locations, components)
    sites = getsites(locations, components)
    allparams = findparams(locations, components, sites)
    finaldict = packageparams(allparams)
    return finaldict
end

function findparams(locations::Array{String,1}, components::Array{String,1}, sites::Array{Array{String,1},1})
    # Returns dictionary with all parameters in their respective arrays.
    checkfor_clashingheaders(locations)
    allparams = Dict()
    for location in reverse(locations)
        type = readtype(location)
        headerparams = readheader(location)
        for headerparam in headerparams
            !haskey(allparams, headerparam) &&
                allparams[headerparam] = createemptyparamsarray(type, components, sites)
            matches = searchfor(type, location, headerparam, components, sites)
            if type == likedata
                for (component, value) in matches
                    idx = findfirst(isequal(component), components)
                    if allparams[headerparam][idx] == 0
                        allparams[headerparam][idx] = value
                    end
                end
            end
            if type == pairdata
                if allparams[headerparams] <: Array{<:Any,1}
                    allparams[headerparams] = convertliketopair(allparams[headerparams])
                for (componentpair, value) in matches
                    idx1 = findfirst(isequal(componentpair[1]), components)
                    idx2 = findfirst(isequal(companentpair[2]), components)
                    if allparams[headerparam][idx1,idx2] == 0
                        allparams[headerparam][idx1,idx2] = value
                    end
                end
            end
            if type == assocdata
                for (assocpair, value) in matches
                    idx1 = findfirst(isequal(assocpair[1][1]), components)
                    idx2 = findfirst(isequal(assocpair[1][2]), components)
                    idx21 = findfirst(isequal(assocpair[2][1]), sites)
                    idx22 = findfirst(isequal(assocpair[2][2]), sites)
                    if allparams[headerparam][idx1,idx2][idx21,idx22] == 0
                        allparams[headerparam][idx1,idx2][idx21,idx22] = value
                    end
                end
            end
        end
    end
    return allparams
end

function searchfor(type::ParameterType, location::String, headerparam::String, components::Array{String,1}, sites::Array{Array{String,1},1})
    # Returns dictionary with all matches in a particular file for one parameter.

end

function readtype(location::String)
end

function readheader(location::String)
end

function retrievesources(locations::Array{String,1}, components::Array{String,1})
end

function findsites(locations::Array{String,1}, components::Array{String,1})
end

function createemptyparamsarray(type::ParameterType, components::Array{String,1}, sites::Array{Array{String,1},1})
    # Creates an empty array of the appropriate size.
    componentslength = length(components)
    type == likedata && return zeros(Double, componentslength)
    type == pairdata && return zeros(Double, componentslength, componentslength)
    if type == assocdata
        output = Array{Array{Float64,2},2}(undef, componentslength, componentslength)
        for i in 1:componentslength, j in 1:componentslength
            output[i,j] = zeros(length(sites[i]), length(sites[j]))
        end
        return output
    end
end


function convertliketopair(params::Array{T,1}) where T
    # Returns a diagonal matrix with the given parameters. 
    paramslength = length(params)
    output = zeros(T, paramslength, paramslength)
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
        if lowervalue != 0 && uppervalue != 0 && lowervalue != uppervalue
            error("Dissimilar non-zero entries exist across diagonal.")
        end
        lowervalue != 0 && (matrix[j,i] = lowervalue)
        uppervalue != 0 && (output[i,j] = uppervalue)
    end
end
