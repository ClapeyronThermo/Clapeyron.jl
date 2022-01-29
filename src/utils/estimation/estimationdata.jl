using CSV, Tables

struct EstimationData
    method::Symbol
    condition_names::Vector{Symbol}
    out_names::Vector{Symbol}
    conditions::Vector{Vector{Union{Float64,Missing}}}
    conditions_error::Vector{Vector{Union{Float64,Missing}}}
    outs::Vector{Vector{Union{Float64,Missing}}}
    outs_error::Vector{Vector{Union{Float64,Missing}}}
end

# Copied from database.jl; methods that are common should be refactored into
# separate files.
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

function extract_dataerror(df::CSV.File, csvheaders::Vector{String}, extract_headers::Vector{String})
    data = Vector{Vector{Union{Float64,Missing}}}(undef, length(extract_headers))
    error = Vector{Vector{Union{Float64,Missing}}}(undef, length(extract_headers))
    for (i, header) in enumerate(extract_headers)
        data[i] = Tables.getcolumn(df, Symbol(header))
        if header * "_error" ∈ csvheaders
            error[i] = Tables.getcolumn(df, Symbol(header * "_error"))
        else
            error[i] = Vector{Missing}(undef, length(df))
        end
    end
    return data, error
end

function get_estimationdata(filepaths::Vector{String})
    estimationdata = Vector{EstimationData}()
    for filepath ∈ filepaths
        method = Symbol(strip(getline(filepath, 2), [',']))
        df = CSV.File(filepath; header=3, pool=0, silencewarnings=true)
        csvheaders = String.(Tables.columnnames(df))
        out_headers = String.(SubString.(filter(x -> startswith(x, "out_") && !endswith(x, "_error"), csvheaders), 5))
        condition_headers = filter(x -> !startswith(x, "out_") && !endswith(x, "_error"), csvheaders)
        conditions, conditions_error = extract_dataerror(df, csvheaders, condition_headers)
        outs, outs_error = extract_dataerror(df, csvheaders, "out_" .* out_headers)
        push!(estimationdata, EstimationData(method, Symbol.(condition_headers), Symbol.(out_headers), conditions, conditions_error, outs, outs_error))
    end
    return estimationdata
end

extract_estimationdata(["saturation_p_rhoL.csv"])
