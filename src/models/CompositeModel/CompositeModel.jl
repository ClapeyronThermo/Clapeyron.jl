#=
struct CompositeModel{𝔽,𝕊} <: EoSModel
    components::Vector{String}
    fluid::𝔽
    solid::𝕊
end
=#

"""
    CompositeModel(components;
    gas = BasicIdeal,
    liquid = RackettLiquid,
    saturation = LeeKeslerSat,
    gas_userlocations = String[],
    liquid_userlocations = String[],
    saturation_userlocations = String[]

Model that holds representations of fluid (and/or solid) that aren't evaluated using the helmholtz energy-based approach used in the rest of the library.

It contains a "fluid" and a "solid" field. there are three available representations for a fluid:
- a helmholtz-based EoS
- Fluid Correlations, consisting in a gas model, a correlation for obtaining the saturation pressure, and a liquid model. both gas and liquid models can optionally be helmholtz models too, but correlations for saturated liquid and vapour are also allowed.
- Activity models, consisting of a liquid activity and a model for the fluid. the fluid model can be a helmholtz-based model, or another `CompositeModel` containing correlations.

When the solid field is specified, some properties (like `volume`) start taking in account the solid phase in their calculations. optionally, there are other models that provide specific correlations for SLE equilibria (like `SolidHfus`)

## Examples:
- Saturation pressure calculated using Correlations:
```julia-repl
#rackett correlation for liquids, DIPPR 101 correlation for the saturation pressure, ideal gas for the vapour volume
julia> model = CompositeModel(["water"],liquid = RackettLiquid,saturation = DIPPR101Sat,gas = BasicIdeal)
Composite Model (Correlation-Based) with 1 component:
 Gas Model: BasicIdeal()
 Liquid Model: RackettLiquid("water")
 Saturation Model: DIPPR101Sat("water")

julia> saturation_pressure(model,373.15)
(101260.56298096628, 1.8234039042643886e-5, 0.030639190960720403)
```

- Bubble Pressure, calculated using fluid correlations and a Raoult solver:
```julia-repl
julia> model = CompositeModel(["octane","heptane"],liquid = RackettLiquid,saturation = DIPPR101Sat,gas = BasicIdeal)
Composite Model (Correlation-Based) with 2 components:
 Gas Model: BasicIdeal()
 Liquid Model: RackettLiquid("octane", "heptane")
 Saturation Model: DIPPR101Sat("octane", "heptane")

julia> bubble_pressure(model,300.15,[0.9,0.1])
(2552.3661540464022, 0.00015684158726046333, 0.9777538974501402, [0.7376170278676232, 0.2623829721323768])
```

- Bubble Pressure, using an Activity Model along with another model for fluid properties:
```julia-repl
#using a helmholtz-based fluid
julia> model = CompositeModel(["octane","heptane"],liquid = UNIFAC,fluid = PR)
Composite Model (γ-ϕ) with 2 components:
 Activity Model: UNIFAC{PR{BasicIdeal, PRAlpha, NoTranslation, vdW1fRule}}("octane", "heptane")
 Fluid Model: PR{BasicIdeal, PRAlpha, NoTranslation, vdW1fRule}("octane", "heptane")

julia> bubble_pressure(model,300.15,[0.9,0.1])
(2694.150594740186, 0.00016898441224336215, 0.9239727973658585, [0.7407077952279438, 0.2592922047720562])

#using a correlation-based fluid
julia> fluidmodel = CompositeModel(["octane","heptane"],liquid = RackettLiquid,saturation = DIPPR101Sat,gas = BasicIdeal); 
model2 = CompositeModel(["octane","heptane"],liquid = UNIFAC, fluid = fluidmodel)
Composite Model (γ-ϕ) with 2 components:
 Activity Model: UNIFAC{PR{BasicIdeal, PRAlpha, NoTranslation, vdW1fRule}}("octane", "heptane")
 Fluid Model: FluidCorrelation{BasicIdeal, RackettLiquid, DIPPR101Sat}("octane", "heptane")

julia> bubble_pressure(model2,300.15,[0.9,0.1])
(2551.6008524130893, 0.00015684158726046333, 0.9780471551726359, [0.7378273929683233, 0.2621726070316766])
```

"""
CompositeModel
"""
    RestrictedEquilibriaModel <: EoSModel

Abstract type of models that implement simplifications over the equality of chemical potentials approach for phase equilibria. subtypes of `RestrictedEquilibriaModel` are the `GammaPhi` (activity + gas), `FluidCorrelation` (for fluid phase change and volume correlations) and `SolidCorrelation` (for solid phase change and solid volume correlations)
"""
abstract type RestrictedEquilibriaModel <: EoSModel end

include("FluidCorrelation.jl")
include("SolidCorrelation.jl")
include("GammaPhi.jl")
include("GenericAncEvaluator.jl")
include("SaturationModel/SaturationModel.jl")
include("LiquidVolumeModel/LiquidVolumeModel.jl")
include("PolExpVapour.jl")
include("SolidModel/SolidHfus.jl")
include("bubble_point.jl")
include("dew_point.jl")

function CompositeModel(components;
    liquid = nothing,
    gas = nothing,
    fluid = nothing,
    solid = nothing,
    saturation = nothing,
    melting = nothing,
    sublimation = nothing,
    gas_userlocations = String[],
    liquid_userlocations = String[],
    fluid_userlocations = String[],
    solid_userlocations = String[],
    saturation_userlocations = String[],
    melting_userlocations = String[],
    sublimation_userlocations = String[],
    verbose = false)

    _components = format_components(components)
    #take care of the solid phase first
    if melting == sublimation == nothing
        init_solid = init_model(solid,components,solid_userlocations,verbose)
    else
        init_solid_phase = init_model(solid,components,solid_userlocations,verbose)
        init_melting = init_model(melting,components,melting_userlocations,verbose)
        init_sublimation = init_model(sublimation,components,sublimation_userlocations,verbose)
        init_solid = SolidCorrelation(_components,init_solid_phase,init_melting,init_sublimation)
    end

    _fluid = init_model(fluid,components,fluid_userlocations,verbose)
    if _fluid isa EoSModel && liquid == gas == saturation == nothing
        #case 1: fluid isa EoSModel. no other model is specified
        if !(_fluid isa ActivityModel)
            init_fluid = _fluid
        else
            error("Activity models only represent the liquid phase. Please specify a fluid phase model.")
        end
    elseif fluid == nothing && !isnothing(liquid) && !isnothing(gas) && !isnothing(saturation)
        #case 2: fluid not specified, V,L,sat specified, use FluidCorrelation struct
        init_gas = init_model(gas,components,gas_userlocations,verbose)
        init_liquid = init_model(liquid,components,liquid_userlocations,verbose)
        init_sat = init_model(saturation,components,saturation_userlocations,verbose)
        init_fluid = FluidCorrelation(_components,init_gas,init_liquid,init_sat)
    elseif _fluid !== nothing && !isnothing(liquid) && (gas == saturation == nothing)
        #case 3: liquid activity and a model for the fluid.
        init_liquid = init_model(liquid,components,liquid_userlocations,verbose)
        if init_liquid isa ActivityModel
            #case 3.a, the fluid itself is a composite model. unwrap the fluid field.
            if _fluid isa CompositeModel
                _fluid = EoSVectorParam(_fluid.fluid,_components)
            else
                _fluid = EoSVectorParam(_fluid,_components)
            end
            init_fluid = GammaPhi(_components,init_liquid,_fluid)
        else
            #case 3.b, one alternative is to leave this as an error.
            init_gas = _fluid
            init_sat = _fluid
            init_fluid = FluidCorrelation(_components,init_gas,init_liquid,init_sat)
        end
    elseif !isnothing(liquid) && (fluid == gas == saturation == nothing)
    #legacy case, maybe we are constructing an activity that has a puremodel
    init_liquid = init_model(liquid,components,liquid_userlocations,verbose)
        if init_liquid isa ActivityModel
            if hasfield(typeof(init_liquid),:puremodel)
                pure = init_liquid.puremodel
            else
                pure = init_puremodel(BasicIdeal(),components,userlocations,verbose)
            end
            init_fluid = GammaPhi(_components,init_liquid,pure)
        else
            throw(ArgumentError("Invalid specification for CompositeModel"))
        end
    else
        throw(ArgumentError("Invalid specification for CompositeModel"))
    end
    return CompositeModel(_components,init_fluid,init_solid)
end

function Base.show(io::IO,mime::MIME"text/plain",model::CompositeModel)
    fluid = model.fluid
    solid = model.solid

    print(io,"Composite Model")
    if fluid isa GammaPhi && solid == nothing
        print(io," (γ-ϕ)")
    elseif fluid isa FluidCorrelation && solid == nothing
        print(io," (Correlation-Based)")
    end
    length(model) == 1 && print(io, " with 1 component:")
    length(model) > 1 && print(io, " with ", length(model), " components:")
    if solid !== nothing
        if solid isa SolidCorrelation
            solid.phase !== nothing && print(io,'\n'," Solid Phase Model: ",model.solid)
            solid.melting !== nothing && print(io,'\n'," Melting Model: ",model.melting)
            solid.sublimation !== nothing && print(io,'\n'," Sublimation Model: ",model.saturation)
        else
            print(io,'\n'," Solid Model: ",solid)
        end
    end

    if fluid !== nothing
        if fluid isa GammaPhi
            print(io,'\n'," Activity Model: ",fluid.activity)
            print(io,'\n'," Fluid Model: ",fluid.fluid.model) #on gamma-phi, fluid is an EoSVectorParam
        elseif fluid isa FluidCorrelation
            fluid.gas !== nothing && print(io,'\n'," Gas Model: ",fluid.gas)
            fluid.liquid !== nothing && print(io,'\n'," Liquid Model: ",fluid.liquid)
            fluid.saturation !== nothing && print(io,'\n'," Saturation Model: ",fluid.saturation)
        else
            fluid !== nothing && print(io,'\n'," Fluid Model: ",fluid)
        end
    end
end

"""
    __gas_model(model::EoSModel)

internal function.
provides the model used to calculate gas properties.
normally, this is the identity, but `CompositeModel` has a gas model by itself.
"""
__gas_model(model::EoSModel) = model
__gas_model(model::CompositeModel) = model.fluid
fluid_model(model::CompositeModel) = model.fluid
solid_model(model::CompositeModel) = model.solid

function volume_impl(model::CompositeModel,p,T,z,phase=:unknown,threaded=false,vol0 = nothing)
    if is_liquid(phase) || is_vapour(phase)
        return volume_impl(model.fluid,p,T,z,phase,threaded,vol0)
    elseif is_solid(phase)
        if !isnothing(model.solid)
            return volume_impl(model.solid,p,T,z,phase,threaded,vol0)
        else
            _0 = zero(p+T+first(z)+one(eltype(model)))
            nan = _0/_0
            return nan
        end
    else #phase = :unknown
        #there is a helmholtz energy model in fluid and solid phases.
        #this requires checking evaluating all volumes and checking
        #what value is the correct one via gibbs energies.
        if !(model.fluid isa GammaPhi) && !(model.fluid isa FluidCorrelation) && !(model.solid isa SolidCorrelation)
            return _volume_impl(model,p,T,z,phase,threaded,vol0)
        else
            #TODO: implement these when we have an actual sublimation-melting empiric model.
            throw(error("automatic phase detection not implemented for $(typeof(model))"))
        end
    end
end

function activity_coefficient(model::CompositeModel,p,T,z=SA[1.]; phase = :unknown, threaded=true)
    return activity_coefficient(model.fluid,p,T,z;phase,threaded)
end

function init_preferred_method(method::typeof(saturation_pressure),model::CompositeModel,kwargs)
    return init_preferred_method(saturation_pressure,model.fluid,kwargs)
end

function init_preferred_method(method::typeof(saturation_temperature),model::CompositeModel,kwargs)
    return init_preferred_method(saturation_temperature,model.fluid,kwargs)
end

function saturation_pressure(model::CompositeModel,T,method::SaturationMethod)
    return saturation_pressure(model.fluid,T,method)
end

crit_pure(model::CompositeModel) = crit_pure(model.fluid)

x0_sat_pure(model::CompositeModel,T) = x0_sat_pure(model.fluid)

x0_psat(model::CompositeModel,T) = x0_psat(model.fluid,T)

function saturation_temperature(model::CompositeModel,p,method::SaturationMethod)
    return saturation_temperature(model.fluid,p,method)
end

#defer bubbledew eq to the fluid field

function init_preferred_method(method::typeof(bubble_pressure),model::CompositeModel,kwargs)
    init_preferred_method(method,model.fluid,kwargs)
end

function init_preferred_method(method::typeof(bubble_temperature),model::CompositeModel,kwargs)
    init_preferred_method(method,model.fluid,kwargs)
end

function init_preferred_method(method::typeof(dew_pressure),model::CompositeModel,kwargs)
    init_preferred_method(method,model.fluid,kwargs)
end

function init_preferred_method(method::typeof(dew_temperature),model::CompositeModel,kwargs)
    init_preferred_method(method,model.fluid,kwargs)
end

function bubble_pressure(model::CompositeModel, T, x, method::BubblePointMethod)
    if !(method isa ActivityBubblePressure) && !(model.fluid isa RestrictedEquilibriaModel)
        throw(ArgumentError("$method not supported by $(typeof(model.fluid))"))
    end
    return bubble_pressure(model.fluid, T, x, method)
end

function bubble_temperature(model::CompositeModel, T, x, method::BubblePointMethod)
    if !(method isa ActivityBubbleTemperature) && !(model.fluid isa RestrictedEquilibriaModel)
        throw(ArgumentError("$method not supported by $(typeof(model.fluid))"))
    end
    return bubble_temperature(model.fluid, T, x, method)
end

function dew_pressure(model::CompositeModel, T, x, method::DewPointMethod)
    if !(method isa ActivityDewPressure)  && !(model.fluid isa RestrictedEquilibriaModel)
        throw(ArgumentError("$method not supported by $(typeof(model.fluid))"))
    end
    return dew_pressure(model.fluid, T, x, method)
end

function dew_temperature(model::CompositeModel, T, x, method::DewPointMethod)
    if !(method isa ActivityDewTemperature)  && !(model.fluid isa RestrictedEquilibriaModel)
        throw(ArgumentError("$method not supported by $(typeof(model.fluid))"))
    end
    return dew_temperature(model.fluid, T, x, method)
end

#Michelsen TPFlash and rachford rice tpflash support
function init_preferred_method(method::typeof(tp_flash),model::CompositeModel{<:Any,Nothing},kwargs)
    init_preferred_method(method,model.fluid,kwargs)
end

__tpflash_cache_model(model::CompositeModel{<:Any,Nothing},p,T,z,equilibrium) = __tpflash_cache_model(model.fluid,p,T,z,equilibrium)

function gibbs_solvation(model::CompositeModel,T)
    binary_component_check(gibbs_solvation,model)
    return gibbs_solvation(model.fluid,T)
end


export CompositeModel
