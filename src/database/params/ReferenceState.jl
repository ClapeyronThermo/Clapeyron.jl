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
        if std_type in (:ashrae,:nbp,:iir,:custom,:no_set)
            return new(components,a0,a1,T0,P0,H0,S0,z0,std_type)
        else
            throw(error("invalid specification for ReferenceState."))
        end
    end
end

function ReferenceState(symbol = :no_set,T0 = NaN,P0 = NaN,H0 = NaN,S0 = NaN,z0 = Float64[])
    _H0 = isnan(H0) ? Float64[] : [H0]
    _S0 = isnan(S0) ? Float64[] : [S0]
    ReferenceState(String[],Float64[],Float64[],T0,P0,_H0,_S0,z0,symbol)
end
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
    if _ref.std_type == :no_set
        reference_state_eval(nothing,V,T,z)
    else
        reference_state_eval(_ref,V,T,z)
    end
end

reference_state_eval(ref::Nothing,V,T,z) = zero(T+first(z)+oneunit(eltype(model)))

function reference_state_eval(ref::ReferenceState,V,T,z)
    ā0 = dot(ref.a0,z)
    ā1 = dot(ref.a1,z)
    return (ā0/T + ā1)/sum(z)
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
    ref.std_type == :no_set && return nothing
    
    #=
    z0 is not set. we:
        - split the model
        - set reference states in each pure model
        - agregate the reference states into the main model
    =#


    if all(iszero,ref.z0) 
        pures = split_model(model)
        z0 = FillArrays.Fill(1.0,length(model))
        _set_reference_state!.(pures,z0)
        pure_refs = reference_state.(pures)
        len = length(model)

        comps,T0,P0,H0,S0 = ref.components,ref.T0,ref.P0,ref.H0,ref.S0
        z0 = ref.z0
        if length(z0) == 0
            resize!(z0,len)
            z0 .= 0
        end

        if length(comps) == 0
            resize!(comps,len)
            comps .= model.components
        end

        if length(H0) == 0
            resize!(H0,len)
        elseif length(H0) == 1
            h0 = H0[1]
            resize!(H0,len)
            H0 .= h0
        end

        if length(S0) == 0
            resize!(S0,len)
        elseif length(S0) == 1
            s0 = S0[1]
            resize!(S0,len)
            S0 .= s0
        end

        a0,a1 = ref.a0,ref.a1
        for prop in (a0,a1)
            resize!(prop,len)
        end
        comps .= model.components
        a0 .= only.(getfield.(pure_refs,:a0))
        a1 .= only.(getfield.(pure_refs,:a1))
    else
        _set_reference_state!(model,ref.z0)
    end
    return model
end

function _set_reference_state!(model,z0 = SA[1.0])
    ref = reference_state(model)
    ref === nothing && return nothing
    type = ref.std_type
    type == :no_set && return nothing
    comps = ref.components
    len = length(model)
    resize!(comps,len)
    comps .= model.components
    T0,P0,H0,S0 = ref.T0,ref.P0,ref.H0,ref.S0
    z0 = ref.z0
    a0,a1 = ref.a0,ref.a1
    resize!(a0,len)
    resize!(a1,len)
    resize!(H0,len)
    resize!(S0,len)
    if length(z0) == 0
        resize!(z0,len)
        z0 .= 0
    end
    @show z0
    a0 .= 0
    a1 .= 0
    R = Rgas(model)
    if type == :ashrae
        single_component_check(set_reference_state!,model)
        T_ashrae = 273.15 - 40
        p,vl_ashrae,_ = saturation_pressure(model,T_ashrae)
        S00 = VT_entropy(model,vl_ashrae,T_ashrae,SA[1.0])
        a1 .= S00 ./R
        H00 = VT_enthalpy(model,vl_ashrae,T_ashrae,SA[1.0])
        a0 .= -H00 ./R
        #ASHRAE: h = 0, s = 0 @ -40C saturated liquid
    elseif type == :nbp
        single_component_check(set_reference_state!,model)
        p_nbp = 101325.0
        T_nbp,vl_nbp,_ = saturation_temperature(model,p_nbp)
        S00 = VT_entropy(model,vl_nbp,T_nbp,SA[1.0])
        a1 .= S00/R
        H00 = VT_enthalpy(model,vl_nbp,T_nbp,SA[1.0])
        a0 .= -H00/R
        #NBP: h=0, s=0 for saturated liquid at 1 atmosphere
    elseif type == :iir
        single_component_check(set_reference_state!,model)
        T_iir = 273.15
        p,vl_iir,_ = saturation_pressure(model,T_iir)
        M = molecular_weight(model,SA[1.0]) #kg/mol
        H_iir = 200*M*1000
        S_iir = 1*M*1000
        S00 = VT_entropy(model,vl_iir,T_iir,SA[1.0])
        a1 .= (S00 - S_iir)/R
        H00 = VT_enthalpy(model,vl_iir,T_iir,SA[1.0])
        a0 .= (-H00 + H_iir)/R
        #IIR: h = 200 kJ/kg, s=1 kJ/kg/K at 0C saturated liquid
    elseif type == :custom
        vl = volume(model,p0,T0,z0)
        H_set = first(H0)
        S_set = first(S0)
        S00 = VT_entropy(model,vl,T0,z0) 
        a1 .= (S00 - S_set)/R
        H00 = VT_enthalpy(model,vl,T0,z0)
        a0 .= (-H00 + H_set)/R
    else
        throw(error("invalid specification for ReferenceState."))
    end
end

export ReferenceState
