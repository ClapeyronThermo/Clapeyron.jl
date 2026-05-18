using CSV, Tables

struct EstimationData{𝔽}
    method::𝔽
    species::Vector{String}
    inputs_name::Vector{Symbol}
    outputs_name::Vector{Symbol}
    inputs::Vector{Vector{Union{Float64,Missing}}}
    outputs::Vector{Vector{Union{Float64,Missing}}}
    inputs_error::Vector{Vector{Union{Float64,Missing}}}
    outputs_error::Vector{Vector{Union{Float64,Missing}}}
    inputs_errortype::Vector{Union{Symbol,Nothing}}
    outputs_errortype::Vector{Union{Symbol,Nothing}}
    weights::Vector{Float64}
end

ERRORTYPES = [:error_abs, :error_rel, :error_std]

# Copied from database.jl; methods that are common should be refactored into
# separate files.
function extract_dataerror(df::CSV.File, csvheaders, extract_headers)
    data = Vector{Vector{Union{Float64,Missing}}}(undef, length(extract_headers))
    error = Vector{Vector{Union{Float64,Missing}}}(undef, length(extract_headers))
    errortype = Vector{Union{Symbol,Nothing}}(nothing, length(extract_headers))
    for (i, header) in enumerate(extract_headers)
        data[i] = Tables.getcolumn(df, Symbol(header))
        for errortype_ ∈ ERRORTYPES
            if header * "_" * String(errortype_) ∈ csvheaders
                error[i] = Tables.getcolumn(df, Symbol(header * "_" * String(errortype_)))
                errortype[i] = errortype_
                break
            end
        end
        if isnothing(errortype)
            error[i] = Vector{Missing}(undef, length(df))
        end
    end
    return data, error, errortype
end

"""
    EstimationData
    EstimationData(filepaths)
## Input parameters:
- `filepaths` or `filepaths_weights`: The filepath of the data used in parameter estimation. Optionally, a tuple containing the weights of each dataset.
## Output: 
An `EstimationData` struct with the following fields:
- `method`: The property estimation method which is used to obtain predictions for a given input
- `inputs_name`: The variable names for the inputs
- `outputs_name`: The variable names for the outputs 
- `inputs`: Vector for each input
- `outputs`: Vector for each output
- `weights`: The weight for this particular dataset
## Description
For a given input data set, produce an `EstimationData` struct.
"""
function EstimationData(filepaths)
    filepaths = flattenfilepaths(String[],filepaths)
    estimationdata = Vector{EstimationData}()
    for filepath ∈ filepaths 
        csv_method = read_csv_options(filepath)
        method = getfield(Main,csv_method.estimator)
        species = csv_method.species
        if isempty(species)
            species=["all"]
        end
        df = read_csv(filepath,DefaultOptions,csv_method[:sep])
        csvheaders = String.(Tables.columnnames(df))
        outputs_headers = chop.(String.(filter(x -> startswith(x, "out_") && !any(endswith.(x, "_" .* String.(ERRORTYPES))), csvheaders)), head=4, tail=0)
        inputs_headers = filter(x -> !startswith(x, "out_") && !any(endswith.(x, "_" .* String.(ERRORTYPES))), csvheaders)
        inputs, inputs_error, inputs_errortype = extract_dataerror(
            df, csvheaders, inputs_headers)
        outputs, outputs_error, outputs_errortype = extract_dataerror(
            df, csvheaders, "out_" .* outputs_headers)
        push!(
            estimationdata,
            EstimationData(
                method,
                species,
                Symbol.(inputs_headers),
                Symbol.(outputs_headers),
                inputs,
                outputs,
                inputs_error,
                outputs_error,
                inputs_errortype,
                outputs_errortype,
                [1.]
                )
            )
    end
    return estimationdata
end

function EstimationData(filepaths_weights::Array{Tuple{Float64, String}})
    filepaths = [filepaths_weights[i][2] for i in 1:length(filepaths_weights)]
    weights = [filepaths_weights[i][1] for i in 1:length(filepaths_weights)]

    filepaths = flattenfilepaths(String[],filepaths)
    estimationdata = Vector{EstimationData}()
    for i ∈  1:length(filepaths)
        csv_method = read_csv_options(filepaths[i])
        method = getfield(Main,csv_method.estimator)
        species = csv_method.species
        if isempty(species)
            species=["all"]
        end
        df = read_csv(filepaths[i],DefaultOptions,csv_method[:sep])
        csvheaders = String.(Tables.columnnames(df))
        outputs_headers = chop.(String.(filter(x -> startswith(x, "out_") && !any(endswith.(x, "_" .* String.(ERRORTYPES))), csvheaders)), head=4, tail=0)
        inputs_headers = filter(x -> !startswith(x, "out_") && !any(endswith.(x, "_" .* String.(ERRORTYPES))), csvheaders)
        inputs, inputs_error, inputs_errortype = extract_dataerror(
            df, csvheaders, inputs_headers)
        outputs, outputs_error, outputs_errortype = extract_dataerror(
            df, csvheaders, "out_" .* outputs_headers)
        push!(
            estimationdata,
            EstimationData(
                method,
                species,
                Symbol.(inputs_headers),
                Symbol.(outputs_headers),
                inputs,
                outputs,
                inputs_error,
                outputs_error,
                inputs_errortype,
                outputs_errortype,
                [weights[i]]
                )
            )
    end
    return estimationdata
end

function Base.show(io::IO, mime::MIME"text/plain", data::EstimationData)
    println(io, "EstimationData{:" * String(Symbol(data.method)) * "}:")
    println(io, " Inputs:")
    for input in data.inputs_name
        println(io, "  :" * String(input))
    end
    println(io, " Outputs:")
    firstloop = true
    for output in data.outputs_name
        !firstloop && println(io, "")
        print(io, "  :" * String(output))
        firstloop = false
    end
end

function Base.show(io::IO, data::EstimationData)
    print(io, "EstimationData{:" * String(Symbol(data.method)) * "}")
end
