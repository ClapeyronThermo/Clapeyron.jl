struct ToEstimate
    params::Vector{Symbol}
    indices::Vector{Vector{Tuple{Int64,Int64}}}  # if nothing, use all
    range_indices::Vector{UnitRange{Int64}} #for vector parameters, create corresponding views of the input parameter
    index_type::Vector{Symbol} #store all index types.

    #lower, upper, guess, unflattened. in the case of no indices, we assume all, that is, we recompute those parameters to be of the same size of the model
    lower::Vector{Vector{Float64}}
    upper::Vector{Vector{Float64}}
    guess::Vector{Vector{Float64}}

    factor::Vector{Float64}
    symmetric::Vector{Bool}
    cross_assoc::Vector{Bool}
    recombine::Vector{Bool}
    ignorefield::Vector{Symbol}
    scalar::Base.RefValue{Bool}
end

"""
    ToEstimate
    ToEstimate(params_dict)

## Input parameters:

A dictionary with the following potential entries:
- `params`: The name of the parameter being fitted (`Symbol`)
- `indices`: The index of the parameter being fitted (`Integer` or `Tuple{Integer,Integer}`)
- `factor`: Factor to multiply parameter being fitted to have it in the correct units (`Float64`)
- `symmetric`: For `PairParam`, if the parameter is symmetric or asymmetric (`Bool`)
- `cross_assoc`: For `AssocParam`, if the parameter is for cross-association (`Bool`)
- `recombine`: For `PairParam`, if the combining rules must be applied for unlike interactions (`Bool`)
- `lower`: Lower bound for the parameter (`Float64`)
- `upper`: Upper bound for the parameter (`Float64`)
- `guess`: Initial guess for the parameter (`Float64`)

## Output:

A `ToEstimate` struct

    ## Description

Turns the input parameter dictionary into a `ToEstimate` struct to be used within the parameter estimation.
"""
function ToEstimate(params_dict::AbstractVector{<:Dict{Symbol}})
    params = Symbol[]
    indices = Vector{Vector{Tuple{Int64,Int64}}}(undef,0)
    range_indices = Vector{UnitRange{Int64}}(undef,0)
    index_type = Symbol[]

    lower = Vector{Float64}[]
    upper = Vector{Float64}[]
    guess = Vector{Float64}[]

    factor = Float64[]
    sym = Bool[]
    recombine = Bool[]
    cross_assoc = Bool[]
    ignorefield = Symbol[]
    for dict in params_dict
        push!(params,dict[:param])
        id = get(dict, :indices, :no_specified)

        if id isa AbstractString
            id = Symbol(String(id))
        end

        if id isa AbstractVector && eltype(id) <: Integer
            id = [(Int(i),Int(i)) for i in id]
        end

        if id isa Symbol
            index_type_i = if id == :no_specified
                :no_specified
            elseif id in (:pair,:unlike)
                :unlike
            elseif id in (:pure,:pures,:diagonal,:like)
                :pures
            elseif id == :all
                :all
            else
                error("ToEstimate: cannot determine index type")
            end
        else
            index_type_i = :specified
        end

        push!(index_type,index_type_i)

         #if id is 0,0, we iterate over all parameters.

        if id isa Integer
            push!(indices, [(convert(Int,id),0)])
        elseif id isa Symbol
            push!(indices,[(0,0)])
        elseif id isa NTuple{2,Int}
            push!(indices, [id])
        else
            push!(indices, id)
        end

        push!(range_indices,1:0)
        n = length(last(indices))
        #numeric optional fields
        append!(factor,get(dict, :factor, 1.))

        lower_ = get(dict, :lower, Float64[])
        lower_ isa AbstractVector && length(lower_) > 1 &&  id == (0,0) && throw(error("ToEstimate: cannot use vector lower bound with unspecified indices"))
        lower_ isa Number ? push!(lower,[lower_]) : push!(lower,lower_)

        upper_ = get(dict, :upper, Float64[])
        upper_ isa AbstractVector && length(upper_) > 1 &&  id == (0,0) && throw(error("ToEstimate: cannot use vector upper bound with unspecified indices"))
        upper_ isa Number ? push!(upper,[upper_]) : push!(upper,upper_)

        guess_ = get(dict, :guess, Float64[])
        guess_ isa AbstractVector && length(guess_) > 1 &&  id == (0,0) && throw(error("ToEstimate: cannot use vector guess with unspecified indices"))
        guess_ isa Number ? push!(guess,[guess_]) : push!(guess,guess_)

        _sym = get(dict,:symmetric,true)
        push!(sym,_sym)
        _recombine = get(dict,:recombine,false)
        push!(recombine,_recombine)
        _cross_assoc = get(dict,:cross_assoc,false)
        push!(cross_assoc,_cross_assoc)
    end

    scalar_indices = all(x -> isone(length(x)),indices)

    return ToEstimate(params, indices, range_indices, index_type, lower, upper, guess, factor, sym, cross_assoc, recombine, ignorefield, Ref(scalar_indices))
end

mutable struct EstimationModel{M} <: EstimationUtils.AbstractEstimationModel{M}
    model::M
    toestimate::ToEstimate
end

set_ignorefield!(pkg_estimate,::Nothing) = nothing
set_ignorefield!(pkg_estimate,ignorefield::Symbol) = push!(pkg_estimate.ignorefield,ignorefield)
set_ignorefield!(pkg_estimate,ignorefield::String) = push!(pkg_estimate.ignorefield,Symbol(ignorefield))

function set_ignorefield!(pkg_estimate,data::Union{Tuple,AbstractArray})
    for i in data
        set_ignorefield!(pkg_estimate,i)
    end
end

function __recursive_getparam(model,sym,ignorefield)
    res = []
    param = getparam(model,sym)
    if !isnothing(param)
        push!(res,param)
    else
        paramnames = fieldnames(typeof(model))
        for paramname ∈ paramnames
            submodel_i = getfield(model,paramname)
            if (submodel_i isa EoSModel) & !(paramname in ignorefield)
                submodel_param = recursive_getparam(submodel_i,sym,ignorefield)
                append!(res,submodel_param)
            end
        end
    end
    return res
end

function recalculate_flatten_estimationmodel!(est_model::EstimationModel)
    #we recalculate the range indices, those are the source of truth between the model and the input indices
    model = est_model.model
    toestimate = est_model.toestimate
    indices = toestimate.indices

    range_indices = toestimate.range_indices
    index_type = toestimate.index_type #stores a vector of symbols indicating each index type
    x0,lb,ub = toestimate.guess,toestimate.lower,toestimate.upper
    paramnames = toestimate.params
    nc = length(model)
    k = 0
    #=
    no_specified: on single comp -> (1,1)

    =#
    for i in eachindex(indices)
        #we check all empty indices
        if index_type[i] != :specified
            #empty model. get all indices, for single component models, this will be (1,1). for other models, it will be a list of indices.
            #check that all available parameters are of the same type and size
            all_params = __recursive_getparam(model,paramnames[i],toestimate.ignorefield)
            param_lengths = map(x -> prod(size(x)),all_params)
            param_types = map(parameterless_type,all_params)
            @assert allequal(param_lengths)
            @assert allequal(param_types)
            nparam = first(param_lengths)

            flag = index_type[i]
            if nc == 1 && nparam == 1 #fast path
                indices[i] = [(1,1)]
            elseif flag == :no_specified && nparam > 1
                @info "EstimationModel: param $(paramnames[i]) has length > 1, but indices were not specified. Defaulting to (1,1) (this could change in future versions of Clapeyron)"
                indices[i] = [(1,1)]
            else
                #get all indices and store those
                id_full = __get_all_indices(first(all_params),flag,toestimate.symmetric[i])
                indices[i] = id_full
            end
        end
        id = indices[i]
        nid = length(id)
        range_indices[i] = (1:nid) .+ k
        k += nid
    end

    est_model.toestimate.scalar[] = all(x -> isone(length(x)),indices)
    return nothing
end

__get_all_indices(param::SingleParameter,flag,sym) = [(i,i) for i in 1:length(param.values)]
function __get_all_indices(param::AssocParam,flag,sym)
        flags = (:no_specified,:pures,:all,:unlike)

    if flag == :pures
        i_pure =  findall(x -> x[1] == x[2],param.values.outer_indices)
        return [(i,0) for i in i_pure]
    elseif flag == :unlike
        i_unlike =  findall(x -> x[1] != x[2],param.values.outer_indices)
        return [(i,0) for i in i_unlike]
    else
        [(i,0) for i in 1:length(param.values.values)]
    end
end


function __get_all_indices(param::PairParameter,flag,sym)
    res = Tuple{Int,Int}[]
    n = LinearAlgebra.checksquare(param.values)

    if flag == :pures
        return [(i,i) for i in 1:n]
    end

    for i = 1:n
        if sym
            if flag == :unlike
                for j = 1:i - 1
                    push!(res,(i,j))
                end
            else
                for j = 1:i
                    push!(res,(i,j))
                end
            end
        else
            for j = 1:n
                flag == :unlike && (i == j) && continue
                push!(res,(i,j))
            end
        end
    end
    return res
end

function EstimationModel(model,pkg_estimate::ToEstimate,ignorefield)
    set_ignorefield!(pkg_estimate,ignorefield)
    est_model = EstimationModel(model,pkg_estimate)
    recalculate_flatten_estimationmodel!(est_model)
    return est_model
end

function EstimationModel(model,toestimate::AbstractVector{<:Dict{Symbol}};ignorefield = nothing)
    return EstimationModel(model,ToEstimate(toestimate),ignorefield)
end

#EstimationUtils API

EstimationUtils.set_eos_parameters!(est_model::EstimationModel,values) = _set_eos_parameters!(est_model.model,est_model,values)

function _set_eos_parameters!(model,est_model::EstimationModel,values)
    estimation = est_model.toestimate
    params = estimation.params
    factor = estimation.factor
    sym = estimation.symmetric
    cross_assoc = estimation.cross_assoc
    idx = estimation.indices
    recombine = estimation.recombine
    ranges = estimation.range_indices
    for (i,ri) in enumerate(ranges)
        id = idx[i]
        value = @view values[ri]
        param = getparam(model,params[i])

        for (k,ik) in enumerate(id)
            __modify_param!(param,ik,value[k],factor[i],recombine[i],sym[i],cross_assoc[i])
        end
    end
    paramnames = fieldnames(typeof(model))
    for paramname ∈ paramnames
        submodel_i = getfield(model,paramname)
        if (submodel_i isa EoSModel) & !(paramname in estimation.ignorefield)
            _set_eos_parameters!(submodel_i,est_model,values)
        end
    end
    recombine!(model)
    return model
end

function EstimationUtils.get_eos_parameters(est_model::EstimationModel)
    model = est_model.model
    T = eltype(model)
    n = sum(length,est_model.toestimate.range_indices)
    values = Vector{T}(undef,n)
    found = Vector{Bool}(undef,n)
    found .= false
    values .= NaN
    return get_eos_parameters!(values,found,model,est_model.toestimate)
end

function get_eos_parameters!(values,found,model,estimation::ToEstimate)
    params = estimation.params
    idx = estimation.indices
    ranges = estimation.range_indices
    factor = estimation.factor
    for (i,paramname) in enumerate(params)
        param = getparam(model,paramname)
        if !isnothing(param)
            id = idx[i]
            ri = ranges[i]
            fi = factor[i]
            found_i = @view found[ri]
            values_i = @view values[ri]
            for (k,ik) in enumerate(id)
                if !found_i[k]
                    vik = __get_param(param,ik)
                    found_i[k] = true
                    values_i[k] = vik/fi
                end
            end
        end
    end

    for paramname ∈ fieldnames(typeof(model))
        submodel_i = getfield(model,paramname)
        if (submodel_i isa EoSModel) & !(paramname in estimation.ignorefield)
            get_eos_parameters!(values,found,submodel_i,estimation)
        end
    end
    if !all(found)
        throw(error("get_eos_parameters: failed to find all parameters."))
    end
    return values
end

function __modify_param!(current_param::Nothing,id::NTuple{2,Int},val,f,recomb,sym,cross_assoc)
    nothing
end


function __modify_param!(current_param::Nothing,id::AbstractVector,val,f,recomb,sym,cross_assoc)
    nothing
end

function __modify_param!(current_param::SingleParameter,id::NTuple{2,Int},val,f,recomb,sym,cross_assoc)
    current_param[id[1]] = val*f
end

function __modify_param!(current_param::PairParameter,id::NTuple{2,Int},val,f,recomb,sym,cross_assoc)
    id1,_id2 = id[1],id[2]
    id2 = iszero(_id2) ? id1 : _id2 #if we use single parameter indices, transform to diagonal indices
    current_param[id1,id2,sym] = val*f
    if (id1==id2) & recomb
        current_param.ismissingvalues[id1,:] .= true
        current_param.ismissingvalues[:,id1] .= true
    elseif id1!=id2
        current_param.ismissingvalues[id1,id2] = false
        current_param.ismissingvalues[id2,id1] = false
    end
end

function __modify_param!(current_param::AssocParam,id::NTuple{2,Int},val,f,recomb,sym,cross_assoc)
    k = id[1]
    if cross_assoc
        ij = current_param.values.outer_indices[k]
        ab = current_param.values.inner_indices[k]
        i,j = ij
        a,b = ab
        ij_mat = current_param.values[i,j]
        ij_mat[a,b] = val*f
        ij_mat[b,a] = val*f
    else
        current_param.values.values[k] = val*f
    end
end

function __modify_param!(current_param::Union{SingleParameter,PairParameter,AssocParam},id::AbstractVector,val,f,recomb,sym,cross_assoc)
    for id_j in id
        modify_param!(current_param,id_j,val,f,recomb,sym,cross_assoc)
    end
end

function __get_param(current_param::SingleParameter,I::NTuple{2,Int})
    current_param[I[1]]
end

function __get_param(current_param::PairParameter,I::NTuple{2,Int})
    i1,_i2 = I[1],I[2]
    i2 = iszero(_i2) ? i1 : _i2
    return current_param[i1,i2]
end

function __get_param(current_param::AssocParam,I::NTuple{2,Int})
    return current_param.values.values[I[1]]
end

function EstimationUtils.symbol_indices(est_model::EstimationModel,syms::Symbol)
    params = est_model.toestimate.params
    sym_ix = findall(isequal(syms),params)
    length(sym_ix) == 0 && throw(error("EstimationUtils.symbol_indices: symbols not found."))
    if length(sym_ix) == 1
        return convert(Vector{Int64},est_model.toestimate.range_indices[sym_ix[1]])
    end
    range_ix = @view est_model.toestimate.range_indices[sym_ix]
    ix = Int[]
    for ri in range_ix
        append!(ix,ri)
    end
    return ix
end

EstimationUtils.parameter_length(est_model::EstimationModel) = sum(length,est_model.toestimate.range_indices)

#=
indexing and broadcasting interface
=#

function __resolve_index(est_model::EstimationModel,ii::Int)
    ranges = est_model.toestimate.range_indices
    indices = est_model.toestimate.indices
    params = est_model.toestimate.params

    if est_model.toestimate.scalar[]
        jk = indices[ii][1]
        return (ii,jk[1],jk[2])
    end

    j = findfirst(x -> ii in x,ranges)
    paramname = params[ii]
    i = isnothing(j) ? 0 : j
    k = ii - first(ranges[i]) + 1

    id_ik = indices[i][k]
    return (i,id_ik[1],id_ik[2])
end

#we suppose that all checks were done before this function
function __get_param(est_model::EstimationModel,ijk::NTuple{3,Int})
    i,j,k = ijk
    __get_param(getfield(est_model.model.params,i),(j,k))
end

__try_get_param(est_model,ii::Int) = __try_get_param(est_model,__resolve_index(est_model,ii))

function __try_get_param(est_model::EstimationModel,ijk::NTuple{3,Int})
    estimation = est_model.toestimate
    model = est_model.model
    i,j,k = ijk
    T = eltype(model)
    i == 0 && (return zero(T)*NaN,false)
    f = estimation.factor[i]
    paramname = estimation.params[i]
    param = getparam(model,paramname)
    if param != nothing
        return __get_param(param,(j,k))/f,true
    end

    for fieldname ∈ fieldnames(typeof(model))
        submodel_i = getfield(model,fieldname)
        if (submodel_i isa EoSModel) & !(paramname in estimation.ignorefield)
            val,success = __try_get_param(model,ijk)
            if success
                return val,success
            end
        end
    end

    return zero(T)*NaN,false
end

function Base.getindex(model::EstimationModel,i::Int)
    val,success = __try_get_param(model,i)
    if success
        return val
    else
        throw(error(lazy"index access to Estimation model failed to get any values."))
    end
end

function Base.setindex!(model::EstimationModel,val,i::Int)
    __modify_param!(model,val,i)
end

__modify_param!(est_model::EstimationModel,val,i::Int) = __modify_param!(est_model.model,est_model,val,__resolve_index(est_model,i))

function __modify_param!(model,est_model::EstimationModel,val,ijk::NTuple{3,Int})
    i,j,k = ijk
    toestimate = est_model.toestimate
    factor = toestimate.factor
    sym = toestimate.symmetric
    cross_assoc = toestimate.cross_assoc
    recombine = toestimate.recombine
    if i != 0
        paramname = toestimate.params[i]
        param = getparam(model,paramname)
        __modify_param!(param,(j,k),val,factor[i],recombine[i],sym[i],cross_assoc[i])
        for paramname ∈ fieldnames(typeof(model))
            submodel_i = getfield(model,paramname)
            if (submodel_i isa EoSModel) & !(paramname in toestimate.ignorefield)
                __modify_param!(submodel_i,est_model,val,ijk)
            end
        end
    end
    return val
end

function __flatten_data(data,toestimate::ToEstimate,default::Float64)
    ranges = toestimate.range_indices
    n = sum(length,ranges)
    flattened_data = fill(default,n)

    for (i,ri) in enumerate(ranges)
        flattened_data_i = @view flattened_data[ri]
        data_i = data[i]
        if length(data_i) == 1
            flattened_data_i .= data_i[1]
        elseif length(data_i) == 0
            #do nothing
        else
            flattened_data_i .= data_i
        end
    end
    return flattened_data
end

EstimationUtils.lower_bounds(model::ToEstimate) = __flatten_data(model.lower,model,-Inf)
EstimationUtils.upper_bounds(model::ToEstimate) = __flatten_data(model.upper,model,Inf)
EstimationUtils.initial_guess(model::ToEstimate) = __flatten_data(model.guess,model,NaN)

function EstimationUtils.initial_guess(est_model::EstimationModel)
    x0 = EstimationUtils.initial_guess(est_model.toestimate)
    all(!isnan,x0) && return x0
    x00 = est_model.toestimate.guess
    ranges = est_model.toestimate.range_indices
    for (i,ri) in enumerate(ranges)
        if length(x00[i]) == 0
            x0i = @view x0[ri]
            param = est_model.toestimate.params[i]
            id_i = est_model.toestimate.indices[i]
            for k in 1:length(ri)
                jk = id_i[k]
                ijk = (i,jk[1],jk[2])
                x0i[k] = __try_get_param(est_model,ijk)[1]
            end
        end
    end
    return x0
end

EstimationUtils.lower_bounds(model::EstimationModel) = EstimationUtils.lower_bounds(model.toestimate)
EstimationUtils.upper_bounds(model::EstimationModel) = EstimationUtils.upper_bounds(model.toestimate)
