"""
    EstimationData{𝔽, L, N, M} <: EstimationUtils.AbstractEstimationLoss

Data structure for parameter estimation.

Fields:
- `method`: function that computes model predictions
- `loss`: loss function comparing prediction and data
- `inputs_name`: names of input variables
- `outputs_name`: names of output variables (from columns prefixed `out_`)
- `inputs`: vector of input tuples (each of length `N`)
- `outputs`: vector of output tuples (each of length `M`)
- `inputs_ismissingvalues`: vector of tuples of length `N` indicating if the corresponding input at the same index is a missing value.
- `outputs_ismissingvalues`:vector of tuples of length `M` indicating if the corresponding output at the same index is a missing value.
- `inputs_error`: vector of error tuples for inputs (NaN where no error)
- `outputs_error`: vector of error tuples for outputs (NaN where no error)
- `inputs_errortype`: error type per input variable (one of `ERRORTYPES` or `:none`)
- `outputs_errortype`: error type per output variable
- `output_weights`: tuple of weights at each evaluation of the `method` function: `sum(w[j]*loss(method[i][j],output[i][j]) for j in 1:M)`. By default is `1.0`
- `data_weights`: vector of weights applied to each data point. by default is 1.0
- `normalize`: if set to `true`, the result of the objective function will be divided by the amount of valid (non-missing) data points.
"""
struct EstimationData{𝔽, L, N, M} <: EstimationUtils.AbstractEstimationLoss
    method::𝔽
    loss::L
    species::Vector{String}
    inputs_name::Vector{Symbol}
    outputs_name::Vector{Symbol}
    inputs::Vector{NTuple{N, Float64}}
    outputs::Vector{NTuple{M, Float64}}
    valid::Vector{Bool}
    inputs_ismissingvalues::Vector{NTuple{N,Bool}}
    outputs_ismissingvalues::Vector{NTuple{M,Bool}}
    inputs_error::Vector{NTuple{N, Float64}}
    outputs_error::Vector{NTuple{M, Float64}}
    inputs_errortype::NTuple{N,Symbol}
    outputs_errortype::NTuple{M,Symbol}
    output_weights::NTuple{M,Float64}
    data_weights::Vector{Float64}
    normalize::Bool
    source::String
end

"""
    EstimationData(file::AbstractString, method = nothing, loss = nothing; output_weights = nothing, input_weights = nothing, species = "all", normalize = true)

Construct an `EstimationData` instance from a Tables.jl‑compatible table, a method and a loss.

- `table_data` must have columns:
    - Input columns with arbitrary names (any column **not** starting with `out_` and not ending with error suffixes).
    - Output columns whose names **start with** `out_`.
    - Optional error columns: for each input/output column `X`, a column named `X_error_abs`, `X_error_rel`, or `X_error_std` (the first found is used).

- `method` is a julia function of the form `f(model,input[1],input[2]...,input[N])::NTuple{M,<:Number}` where N is length of each input tuple and M is the length of each output tuple
- `loss` is a 2-variable loss function (always return a positive value, returns zero when both arguments are equal)
- `output_weights`: tuple of weights at each evaluation of the `method` function: `sum(w[j]*loss(method[i][j],output[i][j]) for j in 1:M)`. By default is `1.0` for each output
- `data_weights`: vector of weights applied to each data point. by default is 1.0 for each data point. They can also be optionally specified in a CSV file.
- `species`: list of species that participate in the evaluation of the objective function. by default, all species are considered.
- `normalize`: if set to `true`, the result of the `objective_function` method will be divided by the amount of valid (non-missing) data points.

If a file argument is passed instead, some properties can be inferred from the CSV options text above the header, in particular:
- The `method` and `loss` keys will be used to get a method and a loss from the `Main` module, if no other method or loss are provided.
- If no `species` argument is provided, the default stored in the CSV `species` key takes priority.
- If no `normalize` argument is provided, the default stored in the CSV `normalize` key takes priority.
- If no `output_weights` argument is provided, the default stored in the CSV `output_weights` key takes priority.
- If no `data_weights` argument is provided, the CSV `data_weights_col` will be used, if available.
"""
function EstimationData(table_data, method, loss; output_weights = nothing,data_weights = nothing,species = nothing,normalize = true, source = nothing)
    cols = Tables.columns(table_data)
    colnames = collect(Symbol, Tables.columnnames(cols))
    nrows = length(Tables.getcolumn(cols, first(colnames)))  # assume at least one column

    # Helper: extract a column as Vector{Float64} (missing → NaN)
    get_col_vals(h) = [ismissing(v) ? NaN : Float64(v) for v in Tables.getcolumn(cols, h)]
    get_missing_vals(h) = [ismissing(v) for v in Tables.getcolumn(cols, h)]


    if data_weights isa Symbol && data_weights != Symbol()
        data_weights_col = data_weights
        ix = findfirst(isequal(data_weights),colnames)
        deleteat!(colnames,ix)
        data_w = get_col_vals(data_weights_col)
    elseif isnothing(data_weights) || data_weights == Symbol()
        data_w = fill(1.0,nrows)
    else
        data_w = convert(Vector{Float64},data_weights)
    end

    out_headers = [h for h in colnames if startswith(string(h), "out_")]
    out_pure    = [h for h in out_headers if !estimationdata_is_error_suffix(h)]
    in_pure     = [h for h in colnames if !startswith(string(h), "out_") && !estimationdata_is_error_suffix(h)]

    inputs_name  = in_pure
    outputs_name = out_pure
    N = length(inputs_name)
    M = length(outputs_name)


    # Build input/output tuples
    input_arrays  = [get_col_vals(h) for h in inputs_name]
    output_arrays = [get_col_vals(h) for h in outputs_name]

    input_ismissing_arrays  = [get_missing_vals(h) for h in inputs_name]
    output_ismissing_arrays = [get_missing_vals(h) for h in outputs_name]

    inputs  = NTuple{N,Float64}[ntuple(k -> Float64(input_arrays[k][i]), Val(N)) for i in 1:nrows]
    outputs = NTuple{M,Float64}[ntuple(k -> Float64(output_arrays[k][i]), Val(M)) for i in 1:nrows]

    inputs_ismissingvalues  = NTuple{N,Bool}[ntuple(k -> input_ismissing_arrays[k][i], Val(N)) for i in 1:nrows]
    outputs_ismissingvalues = NTuple{M,Bool}[ntuple(k -> input_ismissing_arrays[k][i], Val(M)) for i in 1:nrows]

    # Initialize errors with NaN
    inputs_error  = fill(ntuple(i -> NaN,Val{N}()),nrows)
    outputs_error = fill(ntuple(i -> NaN,Val{M}()),nrows)

    inputs_errortype_vec  = fill(:none, N)
    outputs_errortype_vec = fill(:none, M)

    inputs_errortype = estimationdata_fill_errors!(inputs_error, inputs_errortype_vec, inputs_name, Val{N}(), colnames, cols)
    outputs_errortype = estimationdata_fill_errors!(outputs_error, outputs_errortype_vec, outputs_name, Val{M}(), colnames, cols)

    # output weights
    out_w = if output_weights === nothing
        ntuple(i -> 1.0,Val{M}())
    elseif output_weights isa Number
        out_w_only = Float64(output_weights)
        ntuple(i -> out_w_only,Val{M}())
    elseif output_weights isa Tuple || output_weights isa AbstractVector
        n_out_w = length(output_weights)
        if n_out_w == 1 && M > 1
            out_w_only = only(output_weights)
            ntuple(i -> out_w_only,Val{M}())
        elseif n_out_w == M
            ntuple(i -> output_weights[i],Val{M}())
        else
            throw(DimensionMismatch("length of output weights is not equal the amount of outputs"))
        end
    else
        throw(error("could not parse output_weights"))
    end

    valid = [!(any(inputs_ismissingvalues[i]) || any(outputs_ismissingvalues[i])) for i in 1:nrows]
    _species = isnothing(species) ? ["all"] : String.(species)

    _source = isnothing(source) ? "" : ifelse(startswith(source,"Clapeyron Estimator"),"",source)
    return EstimationData{typeof(method), typeof(loss), N, M}(
        method, loss, _species,
        inputs_name, outputs_name,
        inputs, outputs, valid,
        inputs_ismissingvalues, outputs_ismissingvalues,
        inputs_error, outputs_error,
        inputs_errortype, outputs_errortype,
        out_w,
        data_w,
        normalize,
        _source)
end

# Distinguish output/input/error columns
estimationdata_is_error_suffix(s::AbstractString) = any(endswith(s, string("_", e)) for e in ERRORTYPES)
estimationdata_is_error_suffix(s::Symbol) = estimationdata_is_error_suffix(string(s))

function estimationdata_fill_errors!(error_vec, error_type_vec, pure_names, ::Val{len_tuple}, colnames, cols) where len_tuple
    get_col_vals(h,x) = [ismissing(v) ? x : Float64(v) for v in Tables.getcolumn(cols, h)]
    for (j, var_name) in enumerate(pure_names)
        for etype in ERRORTYPES
            err_header = Symbol(string(var_name) * "_" * string(etype))
            if err_header in colnames
                err_col = get_col_vals(err_header,0.0)
                nrows = length(err_col)
                for i in 1:nrows
                    tup = error_vec[i]
                    error_vec[i] = NTuple{len_tuple, Float64}(
                        k == j ? err_col[i] : tup[k] for k in 1:len_tuple
                    )
                end
                error_type_vec[j] = etype
                break  # use first matching error type
            end
        end
    end
    return ntuple(i -> error_type_vec[i],Val{len_tuple}())
end

function Base.show(io::IO, mime::MIME"text/plain", data::EstimationData)
    print(io, typeof(data))

    n = length(data.inputs)
    print(io," with ")
    print(io,n)
    print(io," data point")
    n != 1 && print(io,"s")
    println(io,":")

    print(io, " Method: ")
    println(io,data.method)

    print(io, " Loss: ")
    println(io,data.loss)

    print(io, " Inputs: ")
    if length(data.inputs_name) == 0
        println(io,"none")
    else
        show_pairs(io,data.inputs_name,pair_separator = ", ",quote_string = false)
        println(io)
    end

    print(io, " Outputs: ")
    show_pairs(io,data.outputs_name,pair_separator = ", ",quote_string = false)
    #println(io)

    if data.source != ""
        print(io," Source: ")
        println(io,Base.basename(data.source))
    end
end

function Base.show(io::IO, data::EstimationData)
    print(io, typeof(data))
    print(io,"(")
    n = length(data.input)
    print(io," data point")
    n != 1 && print(io,"s")
    print(io,")")
end

function EstimationUtils.objective_function(data::EstimationData{<:Any, <:Any, N, M}, model) where {N, M}
    # Initialise the objective value (type depends on model, e.g. Float64)
    F = zero(eltype(model))

    # Select the appropriate model for the given species
    if data.species == ["all"]
        model_r = model
    else
        idx = comps_in_equilibria(component_list(model), data.species)
        idx .= .!idx
        model_r = each_split_model(model, idx)
    end

    inputs  = data.inputs
    targets = data.outputs
    valid = data.valid
    dweights = data.data_weights                # per‑point weight
    ow = data.output_weights                    # NTuple{M, Float64}

    npoints = length(data.inputs)
    loss = data.loss
    method = data.method
    #valid_eval = false
    for (inp, targ, dw, valid_point) in zip(inputs, targets, dweights, valid)
        # Compute prediction tuple of length M
        if valid_point
            #valid_eval = true
            pred = method(model_r, inp...)
            # Accumulate weighted loss over the M output variables
            point_loss = sum(ow[j] * loss(pred[j], targ[j]) for j in 1:M)
            F += dw * point_loss
        end
        !isfinite(F) && break
    end

    # normalize if required
    if data.normalize
        F = F / npoints
    end


    #=
    #TODO: do something if no point is valid
    if !valid_eval

    end
    =#

    return F
end

#=

CSV Reading

=#

function read_estimationdata_options(filepath::AbstractString)
    return _read_estimationdata_options(getline(filepath, 2))
end

function _read_estimationdata_options(line::String)
    vec_re = r"\[.*\]"
    maybe_opts_vec = match(vec_re,line)
    json_re = r"\{.*\}"
    maybe_opts_json = match(json_re,line)

    has_csv_options_vec = !isnothing(maybe_opts_vec)
    has_csv_options_json = !isnothing(maybe_opts_json)
    if has_csv_options_json
        __get_estimationdata_options(maybe_opts_json.match,:json)
    elseif has_csv_options_vec
        opts = chop(maybe_opts_vec.match,head = 1,tail = 1)
        return __get_estimationdata_options(opts,:vec)
    else
        return ESTIMATIONDATA_CSV_OPTIONS
    end
end

const ESTIMATIONDATA_CSV_OPTIONS = (loss = Symbol(), method = Symbol(),output_weights = [1.0],data_weights_col = Symbol(),sep = :comma, species = ["all"], normalize = true)

function __get_estimationdata_options(data,type)
    a,b,c,d,e,f,g = "loss","method","output_weights","data_weights_col","sep","species","normalize"
    if type == :vec
        opts_dict = parse_bracket_format(data,true)
        _loss = if haskey(opts_dict,a)
            Symbol(opts_dict[a][2])
        else
            Symbol()
        end

        _method = if haskey(opts_dict,b)
            Symbol(opts_dict[b][2])
        else
            Symbol()
        end

        _output_weights = if haskey(opts_dict,c)
            is_vec,val = opts_dict[c]
            if (!is_vec && count(isspace,val) > 1) || is_vec
                parse.(Float64,_parse_vec(val))
            else
                [parse(Float64,val)]
            end
        else
            [1.0]
        end

         _data_weights_col = if haskey(opts_dict,d)
            Symbol(opts_dict[d][2])
        else
            Symbol()
        end

        _sep = if haskey(opts_dict,e)
            Symbol(opts_dict[e][2])
        else
            :comma
        end

        _species = if haskey(opts_dict,f)
            is_vec,val = opts_dict[f]
            if (!is_vec && count(isspace,val) > 1) || is_vec
                convert(Vector{String},_parse_vec(val))
            else
                [String(val)]
            end
        else
            ["all"]
        end

        _normalize = if haskey(opts_dict,g)
            Base.parse_bool_env(g,String(opts_dict[g][2]),throw = true)
        else
            true
        end

        return (loss = _loss, method = _method,output_weights = _output_weights,data_weights_col = _data_weights_col,sep = _sep, species = _species, normalize = _normalize)
    elseif type == :json
        json_dict = JSON.parse(data)
        _loss = if haskey(json_dict,a)
            Symbol(json_dict[a])
        else
            Symbol()
        end

        _method = if haskey(json_dict,b)
            Symbol(json_dict[b])
        else
            Symbol()
        end

        _output_weights = if haskey(json_dict,c)
            val = json_dict[c]
            if val isa Number
                [Float64(val)]
            else
                convert(Vector{Float64},val)
            end
        else
            [1.0]
        end

         _data_weights_col = if haskey(json_dict,d)
            Symbol(json_dict[d])
        else
            Symbol()
        end

        _sep = if haskey(json_dict,e)
            Symbol(json_dict[e])
        else
            :comma
        end

        _species = if haskey(json_dict,f)
            val = json_dict[f]
            if val isa AbstractString
                [String(val)]
            else
                convert(Vector{String},val)

            end
        else
            ["all"]
        end

        _normalize = if haskey(json_dict,g)
            Base.parse_bool_env(g,String(json_dict[g]),throw = true)
        else
            true
        end

        return (loss = _loss, method = _method,output_weights = _output_weights,data_weights_col = _data_weights_col,sep = _sep, species = _species, normalize = _normalize)
    else
        throw(error("Clapeyron.__get_options: invalid type. expected :json or :vec, got $type"))
    end
end

function EstimationData(data::AbstractString,method = nothing, loss = nothing)
    return EstimationData((NaN,data),method,loss)
end

function EstimationData(data::AbstractVector{<:AbstractString},method = nothing, loss = nothing)
    return EstimationData(only(data),method,loss)
end

function EstimationData(data::AbstractVector{Tuple{<: Any,<:AbstractString}},method = nothing, loss = nothing)
    return EstimationData(only(data),method,loss)
end

function EstimationData(data::Tuple{<: Any,<:AbstractString},method = nothing, loss = nothing)
    output_weights,filepath0 = data
    filepaths = flattenfilepaths(String[],filepath0)
    @assert length(filepaths) == 1
    filepath = filepaths[1]

    estimationdata_options = read_estimationdata_options(filepath)
    table_data = read_csv(filepath,DefaultOptions,estimationdata_options.sep)

    _method = if isnothing(method)
        est = estimationdata_options.method
        if est != Symbol()
            getfield(Main,estimationdata_options.method)
        else
            error("no method specified in CSV and no method passed as argument")
        end
    else
        method
    end

    _loss = if isnothing(loss)
        est = estimationdata_options.loss
        if est != Symbol()
            getfield(Main,estimationdata_options.loss)
        else
            __mse
        end
    else
        loss
    end


    _output_weights = if isnan(output_weights)
        estimationdata_options.output_weights
    else
        output_weights
    end

    data_weights = estimationdata_options.data_weights_col

    species = estimationdata_options.species

    normalize = estimationdata_options.normalize

    EstimationData(table_data, _method, _loss; output_weights = _output_weights, data_weights = data_weights , species = species, normalize = normalize, source = filepath)
end

function estimation_data_from_csvs(filepaths::AbstractVector{<:AbstractString}, method = nothing, loss = nothing)
    filepaths = flattenfilepaths(String[],filepaths)
    res = EstimationData[]
    for filepath in filepaths
        push!(res,EstimationData(filepath,method,loss))
    end
    return res
end

function estimation_data_from_csvs(filepaths_weights::AbstractVector{T}, method = nothing, loss = nothing) where T <: Tuple{<:Real,<:AbstractString}
    filepaths = String[]
    weight = Float64[]
    for filepath_weight in filepaths_weights
        filepaths_i = flattenfilepaths(String[],[filepath_weight[2]])
        n = length(filepaths_i)
        append!(filepaths,filepaths_i)
        append!(weight,fill(convert(Float64,filepath_weight[1]),n))
    end
    res = EstimationData[]
    std_filepath_weigths = collect(zip(weight,filepaths))

    for std_filepath_weight in std_filepath_weigths
        push!(res,EstimationData(std_filepath_weight,method,loss))
    end
    return res
end

