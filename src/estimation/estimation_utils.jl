struct EstimationSolution{L,X,M}
    loss::L
    x::X
    model::M
end


"""
    EstimationUtils

Module storing the interface (abstract types and generic functions) for the
`Claepyron.jl` parameter-estimation framework.

Concrete implementations are expected to subtype `AbstractEstimationModel` and
`AbstractEstimationLoss` and to extend the functions declared here.
"""

module EstimationUtils
    using Roots: @set
#mandatory API for AbstractEstimationModel

"""
    AbstractEstimationModel{M}

Abstract supertype for all parameter-estimation model wrappers in Clapeyron.jl.

An `AbstractEstimationModel{M}` wraps a concrete equation-of-state (EoS) model
of type `M` together with the metadata needed to map a flat parameter vector
`Θ` onto that model's fields and back.

# Required interface

Every concrete subtype must implement:

- [`set_eos_parameters!(estimation_model, Θ)`](@ref): write the flat parameter  vector `Θ` into the wrapped EoS model in-place.
- [`get_eos_parameters(estimation_model)`](@ref): extract a flat parameter  vector `Θ` from the wrapped EoS model.
- [`set_model(estimation_model, new_model)`](@ref): return a new instance of  the estimation model that wraps `new_model` instead of the current one.
- [`get_model(estimation_model)`](@ref): return the EoS model currently stored inside the estimation wrapper.

`set_model` and `get_model` have default implementations that assume that the model is stored in a `model` field.

# Optional interface

- [`parameter_length(estimation_model)`](@ref): return the amount of parameters that are being manipulated.
- [`lower_bounds(estimation_model)`](@ref): return a flat vector of lower bounds for `Θ`, with `-Inf` for unconstrained parameters.
- [`upper_bounds(estimation_model)`](@ref): return a flat vector of upper bounds for `Θ`, with `+Inf` for unconstrained parameters.
- [`initial_guess(estimation_model)`](@ref): return a flat vector of starting values for the optimiser, they may be different that the current stored parameters in the EoS model.

-
# Indexing and broadcasting

The abstract type provides fallback implementations of
`Base.getindex`, `Base.setindex!`, and `Base.copyto!` so that an
`AbstractEstimationModel` can be used like a vector of parameter values.
Indexing by `Symbol` or `Vector{Symbol}` is supported when
`symbol_indices` is defined.

# Concrete subtypes

[`Clapeyron.EstimationModel`](@ref) is the standard concrete subtype.
"""
abstract type AbstractEstimationModel{M} end

"""
    set_eos_parameters!(estimation_model::AbstractEstimationModel{M}, Θ)
    set_eos_parameters!(eos_model::M, estimation_model::AbstractEstimationModel{M}, Θ)

Write the flat parameter vector `Θ` into the EoS model stored in `estimation_model`, updating it in-place.

The two-argument form operates on `estimation_model`'s own internal model.
The three-argument form first calls [`set_model`](@ref) to swap in `eos_model` and then delegates to the two-argument form, allowing a foreign EoS instance to be updated using the same index/factor metadata.
"""
function set_eos_parameters! end


"""
    Θ = get_eos_parameters(estimation_model::AbstractEstimationModel{M})
    Θ = get_eos_parameters(eos_model::M, estimation_model::AbstractEstimationModel{M})

Extract a flat `Vector` of parameter values `Θ` from the EoS model stored in
`estimation_model`.

The one-argument form operates on `estimation_model`'s own internal model.
The two-argument form first calls [`set_model`](@ref) to substitute `eos_model` and then extracts parameters using the same index/factor metadata, allowing inspection of a foreign EoS instance without modifying it permanently.
"""
function get_eos_parameters end

"""
    parameter_length(estimation_model::AbstractEstimationModel{M})::Int

returns the length of the parameter vector `Θ`. Defaults to the length of `get_eos_parameters(estimation_model)`
"""
parameter_length(model::AbstractEstimationModel) = length(get_eos_parameters(model))

"""
    symbol_indices(estimation_model::AbstractEstimationModel{M}, syms::Symbol) -> AbstractVector{Int}

Return the integer indices into the flat parameter vector `Θ` that correspond
to the parameter names in `syms` (a `Symbol` or `AbstractVector{Symbol}`).

Implementing this function is optional; it enables named indexing via
`estimation_model[:param_name]`.
"""
function symbol_indices end

function symbol_indices(model,syms::AbstractVector{Symbol})
    mapreduce(Base.Fix1(symbol_indices,model),vcat,syms)
end
"""
    set_model(estimation_model::AbstractEstimationModel{M}, new_model::M) -> AbstractEstimationModel{M}

Return a new estimation model that is identical to `estimation_model` except
that its wrapped EoS model is replaced by `new_model`.

This is used internally to perform parameter extraction or injection on a
different EoS instance while reusing the same index/factor metadata.

The function defaults to setting the `model` field to the new model.
"""
function set_model(est_model,model)
    @set est_model.model = model
end

"""
    get_model(estimation_model::AbstractEstimationModel{M}) -> M

Return the EoS model currently wrapped by `estimation_model`.

The default implementation just accesses the `model` field.
"""
get_model(est_model) = est_model.model

#indexing interface for #AbstractEstimationModel

function get_eos_parameters(other_model::M,model::AbstractEstimationModel{M}) where M
    new_model = set_model(model,other_model)
    get_eos_parameters(new_model)
end

function set_eos_parameters!(other_model::M,model::AbstractEstimationModel{M},Θ) where M
    new_model = set_model(model,other_model)
    set_eos_parameters!(new_model,Θ)
end

#slow fallbacks, but they allow easy definitions
function Base.getindex(model::AbstractEstimationModel,i::Int)
    Θ = get_eos_parameters(model)
    return Θ[i]
end

Base.BroadcastStyle(::Type{<:T}) where T <: AbstractEstimationModel = Broadcast.Style{T}()
Base.size(m::AbstractEstimationModel) = (parameter_length(m),)

function Base.copyto!(model::AbstractEstimationModel,Θ::AbstractVector)
    set_eos_parameters!(model,Θ)
end

function Base.copyto!(model::AbstractEstimationModel,Θ::Base.Broadcast.Broadcasted)
    set_eos_parameters!(model,collect(Θ))
end

function Base.getindex(model::AbstractEstimationModel,s::Symbol)
    I = symbol_indices(model,s)
    Θ = get_eos_parameters(model)
    return Θ[I]
end

function Base.getindex(model::AbstractEstimationModel,s::AbstractVector{Symbol})
    I = symbol_indices(model,s)
    Θ = get_eos_parameters(model)
    return Θ[I]
end

function Base.setindex!(model::AbstractEstimationModel,Θ_new,i::Int)
    Θ = get_eos_parameters(model)
    Θ[i] = Θ_new
    set_eos_parameters!(model,Θ)
    return Θ_new
end

function Base.setindex!(model::AbstractEstimationModel,Θ_new,s::Symbol)
    i = symbol_indices(model,s)
    Θ = get_eos_parameters(model)
    Θi =  @view Θ[i]
    Θ[i] .= Θ_new
    set_eos_parameters!(model,Θ)
    return Θ_new
end

function Base.setindex!(model::AbstractEstimationModel,Θ_new,s::AbstractVector{Symbol})
    i = symbol_indices(model,s)
    Θ = get_eos_parameters(model)
    Θi =  @view Θ[i]
    Θ[i] .= Θ_new
    set_eos_parameters!(model,Θ)
    return Θ_new
end


#mandatory API for AbstractEstimationLoss
"""
    AbstractEstimationLoss

Abstract supertype for all loss/data containers used in parameter estimation.

A concrete subtype must implement at least [`objective_function`](@ref)

[`Clapeyron.EstimationData`](@ref) is the standard concrete subtype.
"""
abstract type AbstractEstimationLoss end

"""
    objective_function(loss::AbstractEstimationLoss, model) -> Number

Evaluate the scalar objective (loss) for `model` given the experimental data stored in `loss`.
"""
function objective_function end


"""
    lower_bounds(model) -> Vector{Float64}

Return a flat vector of lower bounds for the parameter vector `Θ` used in the estimation procedure.

Elements are `-Inf` for parameters with no lower bound. Used to configure box-constrained optimisers.
"""
function lower_bounds(model::AbstractEstimationModel)
    n = parameter_length(model)
    T = eltype(get_model(model))
    fill(n,T(-Inf))
end

"""
    upper_bounds(model)

Return a flat vector of upper bounds for the parameter vector `Θ` used in the estimation procedure.

Elements are `+Inf` for parameters with no upper bound. Used to configure box-constrained optimisers.
"""
function upper_bounds(model::AbstractEstimationModel)
    n = parameter_length(model)
    T = eltype(get_model(model))
    fill(n,T(Inf))
end

"""
    initial_guess(model)

Return a flat vector of initial parameter values for the optimiser.
The default implementation just returns `get_eos_parameters(model)`.
"""
initial_guess(model::AbstractEstimationModel) = get_eos_parameters(model)

function Base.show(io::IO,::MIME"text/plain",model::AbstractEstimationModel)
    print(io,typeof(model).name.name)
    print(io," for ")
    print(io,get_model(model))
    np = parameter_length(model)
    print(io," with ")
    print(io,np)
    print(io," parameter")
    np != 1 && print(io,"s")
    println(io,":")
    Base.print_matrix(io,EstimationUtils.get_eos_parameters(model))
end

export AbstractEstimationLoss,AbstractEstimationModel
export set_eos_parameters!, get_eos_parameters, get_model, set_model, symbol_indices
export lower_bounds, upper_bounds, initial_guess
export objective_function

end #module

#default loss
__mse(pred,exp) = abs2((pred-exp)/exp)

#default error types used by EstimationData
const ERRORTYPES = [:error_abs, :error_rel, :error_std]
