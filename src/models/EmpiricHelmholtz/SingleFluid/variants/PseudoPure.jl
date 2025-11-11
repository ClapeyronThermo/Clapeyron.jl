#at the implementation level, pseudo-pure fluid models are different from pure fluid models
#just by storing an additional saturation pressure ancillary and by setting pseudo_pure = true.

struct EmpiricPseudoPure <: EmpiricHelmholtzModel
    model::SingleFluid
end


"""
    EmpiricPseudoPure(components;
            userlocations = String[],
            ancillaries = nothing,
            ancillaries_userlocations = String[],
            estimate_pure = false,
            coolprop_userlocations = true,
            Rgas = nothing,
            verbose = false)

## Input parameters
- JSON data (CoolProp and teqp format)

## Input models
- `ancillaries`: a model that provides initial guesses for saturation calculations. if `nothing`, then they will be parsed from the input JSON.

## Description

Instantiates a pseudo-pure Empiric EoS model. Pseudo-pure Empiric models are multicomponent models at a constant composition, represented by one single-component Multiparameter equation.
The equations are only valid in the single phase region up until the saturation line.
Two-phase properties are not available for `EmpiricPseudoPure` models.

## Examples

```julia-repl
julia> model = EmpiricPseudoPure("R410A")
MultiParameter Pseudo-Pure Equation of state for R410A:
 Polynomial power terms: 5
 Exponential terms: 16
"""
function EmpiricPseudoPure(components;
        userlocations = String[],
        ancillaries = nothing,
        ancillaries_userlocations = String[],
        estimate_pure = false,
        coolprop_userlocations = true,
        Rgas = nothing,
        verbose = false,
        idealmodel = nothing,
        ideal_userlocations = String[])

    allow_pseudo_pure = true
    model = SingleFluid(components;userlocations,ancillaries,ancillaries_userlocations,estimate_pure,coolprop_userlocations,Rgas,verbose,idealmodel,ideal_userlocations,allow_pseudo_pure)
    return EmpiricPseudoPure(model)
end

function Base.show(io::IO,mime::MIME"text/plain",model::EmpiricPseudoPure)
    println(io,"MultiParameter Pseudo-Pure Equation of state for $(model.model.components[1]):")
    show_multiparameter_coeffs(io,model.model.residual)
end

Base.length(::EmpiricPseudoPure) = 1

#zero-arg
for f in (:is_splittable,:SingleFluidIdeal,:recombine_impl!,:mw,:idealmodel,:Rgas,:crit_pure,:has_fast_crit_pure)
    @eval begin
        $f(model::EmpiricPseudoPure) = $f(model.model)
    end
end

#1-arg
for f in (:reduced_a_ideal,:T_scale,:p_scale)
    @eval begin
        $f(model::EmpiricPseudoPure,x1) = $f(model.model,x1)
    end
end

lb_volume(model::EmpiricPseudoPure,T,z) = lb_volume(model.model,T,z)


#v,t,z
for f in (:a_ideal,:a_res,:eos_impl,:x0_volume_liquid_lowT,:x0_volume_liquid)
    @eval begin
        $f(model::EmpiricPseudoPure,V,T,z) = $f(model.model,V,T,z)
    end
end

struct SingleFluidFlashMethod <: FlashMethod end

index_reduction(method::SingleFluidFlashMethod,idx_r) = method

for prop in (:bubble_pressure,:bubble_temperature,:dew_pressure,:dew_temperature)
    @eval begin
        function init_preferred_method(method::typeof($prop),model::EmpiricPseudoPure,kwargs)
            return SingleFluidFlashMethod()
        end
    end
end

function pseudo_pure_bubble_pressure(model,T)
    return _eval_generic_anc(model.model.ancillaries.fluid.saturation.data,T)
end

function pseudo_pure_dew_pressure(model,T)
    return _eval_generic_anc(model.model.ancillaries.fluid.saturation.data2,T)
end

function pseudo_pure_bubble_temperature(model,p)
    return _eval_inverse_generic_anc(model.model.ancillaries.fluid.saturation.data,p)
end

function pseudo_pure_dew_temperature(model,p)
    return _eval_inverse_generic_anc(model.model.ancillaries.fluid.saturation.data2,p)
end


function bubble_pressure_impl(model::EmpiricPseudoPure,T,x,method::SingleFluidFlashMethod)
    pl = pseudo_pure_bubble_pressure(model,T)
    vl = volume(model,pl,T,phase = :l)
    vv = volume(model,pl,T,phase = :v)
    return pl,vl,vv,x
end

function dew_pressure_impl(model::EmpiricPseudoPure,T,y,method::SingleFluidFlashMethod)
    pv = pseudo_pure_dew_pressure(model,T)
    vl = volume(model,pv,T,phase = :l)
    vv = volume(model,pv,T,phase = :v)
    return pv,vl,vv,y
end

function bubble_temperature_impl(model::EmpiricPseudoPure,p,x,method::SingleFluidFlashMethod)
    Tl = pseudo_pure_bubble_temperature(model,p)
    vl = volume(model,p,Tl,phase = :l)
    vv = volume(model,p,Tl,phase = :v)
    return Tl,vl,vv,x
end

function dew_temperature_impl(model::EmpiricPseudoPure,p,y,method::SingleFluidFlashMethod)
    Tv = pseudo_pure_dew_temperature(model,p)
    vl = volume(model,p,Tv,phase = :l)
    vv = volume(model,p,Tv,phase = :v)
    return Tv,vl,vv,y
end

export EmpiricPseudoPure
