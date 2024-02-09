#we take the nomenclature from SingleFluid ideal params:
#a_ideal(V,T) = f(T) + a0 + a1*T
mutable struct ReferenceState <: ClapeyronParam
    components::Vector{String}
    a0::Vector{Float64}
    a1::Vector{Float64}
    T0::Float64
    P0::Float64
    H0::Vector{Float64}
    S0::Vector{Float64}
    z0::Vector{Float64}
    std_type::Symbol
    function ReferenceState(components,a0,a1,T0,P0,H0,S0,z0,std_type)
        if std_type in (:ASHRAE,:NBP,:IIR,:CUSTOM_PURE,:CUSTOM_COMP,:NO_SET)
            return new(components,a0,a1,T0,P0,H0,S0,z0,std_type)
        else
            throw(error("invalid specification for ReferenceState."))
        end
    end
end

ReferenceState() = ReferenceState(String[],Float64[],Float64[],NaN,NaN,Float64[],Float64[],Float64[],:NO_SET)

#by default, the reference state is stored in the idealmodel params. unwrap until
#reaching that
reference_state(model) = reference_state(idealmodel(model))

function reference_state(model::IdealModel)
    return __reference_state(model)
end

@generated function __reference_state(model)
    if hasfield(model,:params)
        params = fieldtype(model,:params)
        if hasfield(params,:reference_state)
            if fieldtype(params,:reference_state) == ReferenceState
                return quote model.params.reference_state end
            end     
        end
    end
    return quote nothing end
end

function reference_state_eval(model::EoSModel,V,T,z)
    _ref = reference_state(model)
    if _ref.std_type == :NO_SET
        reference_state_eval(nothing,V,T,z)
    else
        reference_state_eval(_ref,V,T,z)
    end
end

reference_state_eval(ref::Nothing,V,T,z) = zero(T+first(z)+oneunit(eltype(model)))

function reference_state_eval(ref::ReferenceState,V,T,z)
    ā0 = dot(ref.a0,z)
    ā1 = dot(ref.a0,z)
    return (ā0 + ā1*T)/sum(z)
end

has_reference_state(x) = false
has_reference_state(model::EoSModel) = has_reference_state(typeof(model))
function has_reference_state(::Type{M}) where M <: EoSModel
    return hasfield(M,:params) && has_reference_state(fieldtype(M,:params))
end

function has_reference_state(model::Type{T}) where T
    return hasfield(T,:reference_state) && (fieldtype(T,:reference_state) == ReferenceState)
end

function set_reference_state!(model::EoSModel)
    ref = reference_state(model)
    #handle cases where we don't need to do anything
    ref === nothing && return nothing
    ref.std_type == :NO_SET && return nothing
    if all(iszero,)
    #TODO: actually implement the fit
end

export ReferenceState
