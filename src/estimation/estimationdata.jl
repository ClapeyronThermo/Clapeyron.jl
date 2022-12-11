using CSV, Tables

struct EstimationData{𝔽}
    method::𝔽
    inputs_name::Vector{Symbol}
    outputs_name::Vector{Symbol}
    inputs::Vector{Vector{Union{Float64,Missing}}}
    outputs::Vector{Vector{Union{Float64,Missing}}}
    inputs_error::Vector{Vector{Union{Float64,Missing}}}
    outputs_error::Vector{Vector{Union{Float64,Missing}}}
    inputs_errortype::Vector{Union{Symbol,Nothing}}
    outputs_errortype::Vector{Union{Symbol,Nothing}}
end

ERRORTYPES = [:error_abs, :error_rel, :error_std]

# Copied from database.jl; methods that are common should be refactored into
# separate files.
function extract_dataerror(df::CSV.File, csvheaders::Vector{String}, extract_headers::Vector{String})
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

function EstimationData(filepaths::Vector{String})
    filepaths = flattenfilepaths(String[],filepaths)
    estimationdata = Vector{EstimationData}()
    for filepath ∈ filepaths 
        csv_method = read_csv_options(filepath).estimator
        method = getfield(Main,csv_method)
        df = read_csv(filepath)
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
                Symbol.(inputs_headers),
                Symbol.(outputs_headers),
                inputs,
                outputs,
                inputs_error,
                outputs_error,
                inputs_errortype,
                outputs_errortype
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
