struct ToEstimate{M}
    model::M
    params::Vector{Symbol}
    indices::Vector{Union{Int,Tuple{Int,Int},Vector{Int},Nothing}}  # if nothing, use all
    factor::Vector{Union{Float64,Nothing}}
    lower::Vector{Union{Vector{Union{Float64,Nothing}},Nothing}}
    upper::Vector{Union{Vector{Union{Float64,Nothing}},Nothing}}
    guess::Vector{Union{Vector{Union{Float64,Nothing}},Nothing}}  # if nothing, use current
    symmetric::Vector{Bool}
    cross_assoc::Vector{Bool}
    recombine::Vector{Bool}
    ignorefield::Vector{Symbol}
end

"""
    ToEstimate
    ToEstimate(params_dict)
    ToEstimate(model,params_dict)

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
function ToEstimate(params_dict::Vector{Dict{Symbol,Any}})
    params = Vector{Symbol}(undef,0)
    indices = Vector{Union{Integer,Tuple{Integer,Integer},Vector,Nothing}}(nothing,0)
    factor = Vector{Union{Float64,Nothing}}(nothing,0)
    lower = Float64[]
    upper = Float64[]
    guess = Float64[]
    sym = Vector{Bool}(undef,0)
    recombine = Vector{Bool}(undef,0)
    cross_assoc = Vector{Bool}(undef,0)
    ignorefield = Vector{Symbol}(undef,0)
    for dict in params_dict
        push!(params, dict[:param])
        id = get(dict, :indices, (1,1))
        n = (id isa Tuple{Int,Int} || id isa Int) ? 1 : length(id)
        push!(indices, id)
        push!(factor, get(dict, :factor, 1.))
        lower_ = get(dict, :lower, fill(-Inf,n))
        append!(lower, lower_)
        upper_ = get(dict, :upper, fill(Inf,n))
        append!(upper_, lower_)
        guess_ = get(dict, :guess, fill(NaN,n))
        append!(guess,guess_)
        _sym = get(dict,:symmetric,true)
        push!(sym,_sym)
        _recombine = get(dict,:recombine,false)
        push!(recombine,_recombine)
        _cross_assoc = get(dict,:cross_assoc,false)
        push!(cross_assoc,_cross_assoc)
    end
    return ToEstimate(params, indices, factor, lower, upper, guess, sym, cross_assoc, recombine,ignorefield)
end

mutable struct EstimationModel{M}
    model::M
    toestimate::ToEstimate
end

set_ignorefield!(pkg_estimate,::Nothing) = nothing
set_ignorefield!(pkg_estimate,::Symbol) = push!(pkg_estimate.ignorefield,ignorefield)
set_ignorefield!(pkg_estimate,::String) = push!(pkg_estimate.ignorefield,Symbol(ignorefield))
function set_ignorefield!(pkg_estimate,data::Union{Tuple,AbstractArray})
    for i in data
        set_ignorefield!(pkg_estimate,i)
    end
end

function EstimationModel(model,pkg_estimate::ToEstimate,ignorefield)
    set_ignorefield!(pkg_estimate,ignorefield)
    guess_from_data = pkg_estimate.guess
    est_model = EstimationModel(model,pkg_estimate)
    #update guess
    if any(isnan,guess_from_data)
        guess_from_eos = get_eos_parameters(est_model)
        for (i,x0i) in guess_from_data
            if isnan(x0i)
                guess_from_data[i] = guess_from_eos[i]
            end
        end
    end
end

function EstimationModel(model,toestimate::Vector{Dict{Symbol,Any}};ignorefield = nothing)
    return EstimationModel(model,ToEstimate(toestimate),ignorefield)
end

function return_model(estimation::EstimationModel,model::EoSModel,values)
    T = Base.promote_eltype(model,values)
    return_model!(estimation,promote_model(T,model),values)
end

export return_model

"""
    return_model
    return_model(estimation::Estimation,model,params)
    return_model(model::EstimationModel,params)


## Input parameters:
- `estimation`: The estimator object
- `model`: The model whose parameters we are varying
- `params`: The new parameters which we want to change
## Output:
- `model`: The new model with the updated parameters
## Description
Based on the parameters provided and the estimator, a new model is produced from the input.
"""
function return_model(estimation::Estimation,model::EoSModel,values)
    T = Base.promote_eltype(model,values)
    return_model!(estimation,promote_model(T,model),values)
end

function return_model!(
    estimation::Estimation,
    model::EoSModel,
    values)
    return set_eos_parameters!(model,estimation,values)
end

function set_eos_parameters!(model,estimation::EstimationModel,values)
    params = estimation.params
    factor = estimation.factor
    sym = estimation.symmetric
    cross_assoc = estimation.cross_assoc
    idx = estimation.indices
    recombine = estimation.recombine
    if isdefined(model,:params)
        for (i, param) in enumerate(params)
            if isdefined(model.params,param)
                f = factor[i]
                id = idx[i]
                recomb = recombine[i]
                val_i,sym_i,cross_assoc_i = values[i],sym[i],cross_assoc[i]
                current_param = getfield(model.params, param)
                __modify_param!(current_param,id,val_i,f,recomb,sym_i,cross_assoc_i)
            end
        end 
    end
    paramnames = fieldnames(typeof(model))
    for paramname ∈ paramnames
        submodel_i = getfield(model,paramname)
        if (submodel_i isa EoSModel) & !(paramname in estimation.ignorefield)
            set_eos_parameters!(submodel_i,estimation,values)
        end
    end
    recombine!(model)
end

function get_eos_parameters(est_model::EstimationModel)
    model = est_model.model
    T = eltype(model)
    values = Vector{T}(undef,0)
    return get_eos_parameters!(values,model,est_model.toestimate)
end

function get_eos_parameters!(values,est_model,estimation::ToEstimate)
    model = est_model.model
    estimation = est_model.toestimate
    params = estimation.params
    sym = estimation.symmetric
    idx = estimation.indices
    if isdefined(model,:params)
        for (i, param) in enumerate(params)
            if isdefined(model.params,param)
                id = idx[i]
                sym_i,sym[i]
                current_param = getfield(model.params, param)
                if id isa AbstractVector
                    for ik in id
                        vik = __get_param(current_param,ik,sym_i)
                        push!(values,vik)
                    end
                else
                    push!(values,__get_param(current_param,id,sym_i))
                end
            end
        end 
    end
    for paramname ∈ fieldnames(typeof(model))
        submodel_i = getfield(model,paramname)
        if (submodel_i isa EoSModel) & !(paramname in estimation.ignorefield)
            get_eos_parameters!(values,submodel_i,estimation)
        end
    end
    return values
end

function __modify_param!(current_param::SingleParameter,id::Union{Tuple,Integer},val,f,recomb,sym,cross_assoc)
    current_param[id[1]] = val*f
end

function __modify_param!(current_param::PairParameter,id::Union{Tuple,Integer},val,f,recomb,sym,cross_assoc)
    id1,id2 = id[1],id[2]
    current_param[id1,id2,sym] = val*f
    if (id1==id2) & recomb
        current_param.ismissingvalues[id1,:] .= true
        current_param.ismissingvalues[:,id1] .= true
    elseif id1!=id2
        current_param.ismissingvalues[id1,id2] = false
        current_param.ismissingvalues[id2,id1] = false
    end
end

function __modify_param!(current_param::AssocParam,id::Union{Tuple,Integer},val,f,recomb,sym,cross_assoc)
    id1 = id[1]
    current_param.values.values[id1] = val*f
    if cross_assoc
        current_param.values.values[id1+1] = val*f
    end
end

function __modify_param!(current_param::Union{SingleParameter,PairParameter,AssocParam},id::AbstractVector,val,f,recomb,sym,cross_assoc)
    for id_j in id
        modify_param!(current_param,id_j,val,f,recomb,sym,cross_assoc)
    end
end

#==#

function __get_param(current_param::SingleParameter,I::Union{Tuple,Integer},sym)
    current_param[I]
end

function __get_param(current_param::PairParameter,I::Union{Tuple,Integer},sym)
    i1,i2 = I[1],O[2]
    return current_param[i1,i2,sym]
end

function __get_param(current_param::AssocParam,I::Union{Tuple,Integer},sym)
    return current_param.values.values[I[1]]
end

function __get_param(current_param::Union{SingleParameter,PairParameter,AssocParam},I::AbstractVector,sym)
    [__get_param(current_param,i,sym) for i in I]
end

#=
indexing and broadcasting interface
=#

function __getindex_EstimationModel(model::EstimationModel,ii::Union{Int,Symbol})
    estimation = model.toestimate
    params = estimation.params
    sym = estimation.symmetric
    idx = estimation.indices
    model = est_model.model
    if isdefined(model.model,:params)
        if ii == Int
            i = ii
        else
            j = findfirst(isequal(ii),params)
            if isnothing(j)
                i = 0
            else
                i = j
            end
        end
        if i != 0 && isdefined(model.params,params[i])
            id = idx[i]
            sym_i = sym[i]
            current_param = getfield(model.params, param)
            id isa AbstractVector &&  throw(error(lazy"Base.getindex(model::EstimationModel,i) does not work with abstract vectors fields yet"))
            return __get_param(current_param,id,sym_i),true
        end
    end

    for paramname ∈ fieldnames(typeof(model))
        submodel_i = getfield(model,paramname)
        if (submodel_i isa EoSModel) & !(paramname in estimation.ignorefield)
            val,success = __getindex_EstimationModel(model::EstimationModel,ii::Union{Int,Symbol})
            if success
                return val,success
            end
        end
    end
    T = eltype(model)
    return zero(T)*NaN,false
end

function Base.getindex(model::EstimationModel,i::Union{Int,Symbol})
    val,success = __getindex_EstimationModel(model,i)
    if success
        return val
    else
        throw(error(lazy"index access to Estimation model failed to get any indices."))
    end
end

