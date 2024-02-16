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
    phase::Symbol
    std_type::Symbol
end
"""
    ReferenceState(type::Symbol = :no_set;T0 = NaN;P0 = NaN,H0 = NaN,S0 = NaN,phase = :unknown,z0 = Float64[])

Parameter used to define a reference state for enthalpy and entropy, normally stored in the ideal model. 
when set, it calculates a set of `a0` and `a1` values such as the entropy and enthalpy at a specified point are fixed.

the `type` argument accepts the following standalone options:
- `:no_set`: it returns the current defaults stablished by the equation of state.
- `:ashrae`: h = s = 0 at -40C saturated liquid
- `:iir`: h = 200.0 kJ/kg, s=1.0 kJ/kg/K at 0C saturated liquid
- `:nbp`: h = s = 0 at 1 atm saturated liquid

it also accepts the following options, that require additional specifications:
- `:volume` h = H0, s = S0, at T = T0, v = `volume(model,P0,T0,z0,phase = phase)`
- `:saturation_pressure` h = H0, s = S0, at T = T0, saturated phase (specified by the `phase` argument)
- `:saturation_temperature` h = H0, s = S0, at p = P0, saturated phase (specified by the `phase` argument)

If `z0` is not specified, the reference state calculation will be done for each component separately.

## Examples
```
julia> model = PCSAFT(["water","pentane"],idealmodel = ShomateIdeal,reference_state = ReferenceState(:nbp))
PCSAFT{ShomateIdeal, Float64} with 2 components:
 "water"
 "pentane"
Contains parameters: Mw, segment, sigma, epsilon, epsilon_assoc, bondvol

julia> model2 = PCSAFT(["water","pentane"],idealmodel = ShomateIdeal,reference_state = :nbp) #equivalent
PCSAFT{ShomateIdeal, Float64} with 2 components:
 "water"
 "pentane"
Contains parameters: Mw, segment, sigma, epsilon, epsilon_assoc, bondvol

julia> pure = split_model(model)
2-element Vector{PCSAFT{ShomateIdeal, Float64}}:
 PCSAFT{ShomateIdeal, Float64}("water")
 PCSAFT{ShomateIdeal, Float64}("pentane")

julia> T,vl,_ = saturation_temperature(pure[1],101325.0) #saturated liquid at 1 atm
(373.2706553019503, 2.0512186595412677e-5, 0.03006573003253086)

julia> enthalpy(pure[1],101325.0,T)
-5.477897970382323e-6

julia> entropy(pure[1],101325.0,T)
5.009221069190994e-9
```
"""
ReferenceState

function ReferenceState(symbol = :no_set;T0 = NaN,P0 = NaN,H0 = NaN,S0 = NaN,phase = :unknown, z0 = Float64[])
    _H0 = isnan(H0) ? Float64[] : [H0]
    _S0 = isnan(S0) ? Float64[] : [S0]
    _symbol = if !isnan(T0) & !isnan(P0) & (symbol == :no_set)
        :volume
    else
        symbol
    end
    ReferenceState(String[],Float64[],Float64[],T0,P0,_H0,_S0,z0,phase,_symbol)
end

__init_reference_state_kw(::Nothing) = ReferenceState()
__init_reference_state_kw(s::Symbol) = ReferenceState(s)
__init_reference_state_kw(ref::ReferenceState) = deepcopy(ref)

#by default, the reference state is stored in the idealmodel params. unwrap until
#reaching that
"""
    reference_state(model)::Union{ReferenceState,Nothing}

Returns the reference state of the input model, if available. Returns `nothing` otherwise.

## Examples
```julia-repl
julia> reference_state(PCSAFT("water"))
false

julia> has_reference_state(PCSAFT("water",idealmodel = ReidIdeal))
true

julia> reference_state(PCSAFT("water",idealmodel = MonomerIdeal)) #has reference state, it is not initialized.
ReferenceState(String[], Float64[], Float64[], NaN, NaN, Float64[], Float64[], Float64[], :unknown, :no_set)

julia> reference_state(PCSAFT("water",idealmodel = MonomerIdeal, reference_state = ReferenceState(:nbp))) #has an initialized reference state
ReferenceState(["water"], [33107.133379491206], [17.225988924236503], NaN, NaN, [0.0], [0.0], [0.0], :unknown, :nbp)
```
"""
function reference_state end
reference_state(model::EoSModel) = reference_state(idealmodel(model))
reference_state(::Nothing) = nothing

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
    reference_state_eval(_ref,V,T,z)    
end

reference_state_eval(ref::Nothing,V,T,z) = zero(1.0*T+first(z))

function reference_state_eval(ref::ReferenceState,V,T,z)
    if ref.std_type == :no_set
        return zero(1.0*T + first(z))
    end
    ā0 = dot(ref.a0,z)
    ā1 = dot(ref.a1,z)
    return ā0 + ā1*T
end

"""
    has_reference_state(model)::Bool

Checks if the input `EoSModel` has a reference state. Returns `true/false`

## Examples

```julia-repl
julia> has_reference_state(PCSAFT("water"))
false

julia> has_reference_state(PCSAFT("water",idealmodel = ReidIdeal))
true
```

Note that the default idealmodel (`BasicIdeal`) does not allow for setting reference states.
"""
has_reference_state(x) = !isnothing(reference_state(x))
has_reference_state(x::Type{T}) where T = has_reference_state_type(T)

function has_reference_state_type(::Type{model}) where model
    if hasfield(model,:params)
        params = fieldtype(model,:params)
        if hasfield(params,:reference_state)
            if fieldtype(params,:reference_state) == ReferenceState
                return true
            end     
        elseif hasfield(model,:idealmodel)
            return has_reference_state_type(fieldtype(model,:idealmodel))
        end
    elseif hasfield(model,:idealmodel)
        return has_reference_state_type(fieldtype(model,:idealmodel))
    end
    return false
end

function set_reference_state!(model::EoSModel;verbose = false)
    #handle cases where we don't need to do anything
    ref = reference_state(model)
    ref === nothing && return nothing
    ref.std_type == :no_set && return nothing
    if verbose
        @info "Calculating reference states for $model..."
        @info "Reference state type: $(info_color(ref.std_type))"
    end

    #allocate the appropiate caches.
    initialize_reference_state!(model,ref)
    if all(iszero,ref.z0) #pure case
        pures = split_model(model)
        _set_reference_state!.(pures)
        pure_refs = reference_state.(pures)
        ref.a0 .= only.(getfield.(pure_refs,:a0))
        ref.a1 .= only.(getfield.(pure_refs,:a1))
    else
        _set_reference_state!(model,ref.z0)
    end
    return model
end

function _set_reference_state!(model,z0 = SA[1.0],ref = reference_state(model))
    ref === nothing && return nothing
    type = ref.std_type
    type == :no_set && return nothing
    
    T0,P0,H0,S0 = ref.T0,ref.P0,ref.H0,ref.S0
    a0,a1 = ref.a0,ref.a1
    R = Rgas(model)
    
    if type == :ashrae
        #ASHRAE: h = 0, s = 0 @ -40C saturated liquid
        single_component_check(set_reference_state!,model)
        T_ashrae = 273.15 - 40
        _a0,_a1 = calculate_reference_state_consts(model,:saturation_pressure,T_ashrae,NaN,0.,0.,SA[1.0],:liquid)
        a0 .= _a0
        a1 .= _a1
    elseif type == :nbp
        #NBP: h=0, s=0 for saturated liquid at 1 atmosphere
        single_component_check(set_reference_state!,model)
        p_nbp = 101325.0
        _a0,_a1 = calculate_reference_state_consts(model,:saturation_temperature,NaN,p_nbp,0.,0.,SA[1.0],:liquid)
        a0 .= _a0
        a1 .= _a1
    elseif type == :iir
        #IIR: h = 200 kJ/kg, s=1 kJ/kg/K at 0C saturated liquid
        single_component_check(set_reference_state!,model)
        T_iir = 273.15
        M = molecular_weight(model,SA[1.0]) #kg/mol
        H_iir = 200*M*1000
        S_iir = 1*M*1000
        _a0,_a1 = calculate_reference_state_consts(model,:saturation_pressure,T_iir,NaN,H_iir,S_iir,SA[1.0],:liquid)
        a0 .= _a0
        a1 .= _a1
    elseif type in (:volume,:saturation_pressure,:saturation_temperature)
        _a0,_a1 = calculate_reference_state_consts(model,type,T0,P0,first(H0),first(S0),z0,ref.phase)
        a0 .= _a0
        a1 .= _a1
    else
        throw(error("invalid specification for ReferenceState."))
    end
end

function initialize_reference_state!(model,ref = reference_state(model))
    comps,T0,P0,H0,S0 = ref.components,ref.T0,ref.P0,ref.H0,ref.S0
    z0 = ref.z0
    len = length(model)
    pure_check = length(z0) == 0

    if pure_check
        resize!(z0,len)
        z0 .= 0
    end

    if length(comps) == 0
        resize!(comps,len)
        comps .= model.components
    else
        #this means the ReferenceState struct was already initialized. check for inconsistencies in size
        check_arraysize(model,comps)
    end

    if length(H0) == 0
        resize!(H0,len)
        H0 .= 0
    elseif length(H0) == 1
        h0 = H0[1]
        resize!(H0,len)
        H0 .= h0
    end

    if length(S0) == 0
        resize!(S0,len)
        S0 .= 0
    elseif length(S0) == 1
        s0 = S0[1]
        resize!(S0,len)
        S0 .= s0
    end
    resize!(ref.a0,len)
    resize!(ref.a1,len)

    if !pure_check
        h0 = H0[1]
        if !all(isequal(h0),H0)
            throw(ArgumentError("cannot set enthalpy to different values when evaluating in a multicomponent reference state."))
        end

        s0 = S0[1]
        if !all(isequal(s0),S0)
            throw(ArgumentError("cannot set entropy to different values when evaluating in a multicomponent reference state."))
        end
    end
    return ref
end

function calculate_reference_state_consts(model,type,T0,P0,H0,S0,z0,phase)
    if type == :saturation_pressure
        p,vl,vv = saturation_pressure(model,T0)
        v = is_liquid(phase) ? vl : vv
        T = T0
    elseif type == :saturation_temperature
        T,vl,vv = saturation_temperature(model,P0)
        v = is_liquid(phase) ? vl : vv
    elseif type == :volume
        v = volume(model,P0,T0,z0,phase = phase)
        T = T0
    else
    end
    return __calculate_reference_state_consts(model,v,T,z0,H0,S0)
end

function __calculate_reference_state_consts(model,v,T,z,H0,S0)
    ∑z = sum(z)
    S00 = VT_entropy(model,v,T,z)
    a1 = (S00 - S0)#/∑z
    H00 = VT_enthalpy(model,v,T,z)
    a0 = (-H00 + H0)#/∑z
    return a0,a1
end

export ReferenceState,reference_state,has_reference_state,set_reference_state!
