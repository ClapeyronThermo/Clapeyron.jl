#at the implementation level, pseudo-pure fluid models are different from pure fluid models
#just by storing an additional saturation pressure ancillary and by setting pseudo_pure = true.

struct SingleFluidFlashMethod <: FlashMethod end

index_reduction(method::SingleFluidFlashMethod,idx_r) = method
is_pseudo_pure(model::SingleFluid) = model.properties.pseudo_pure

for prop in (:bubble_pressure,:bubble_temperature,:dew_pressure,:dew_temperature)
    @eval begin
        function init_preferred_method(method::typeof($prop),model::SingleFluid,kwargs)
            return SingleFluidFlashMethod()
        end
    end
end

function bubble_pressure_impl(model,T,x,method::SingleFluidFlashMethod)
    pl = _eval_generic_anc(model.ancillaries.fluid.saturation.data,T)
    vl = volume(model,pl,T,phase = :l)
    vv = volume(model,pl,T,phase = :v)
    return pl,vl,vv,x
end

function dew_pressure_impl(model,T,y,method::SingleFluidFlashMethod)
    pv = _eval_generic_anc(model.ancillaries.fluid.saturation.data2,T)
    vl = volume(model,pv,T,phase = :l)
    vv = volume(model,pv,T,phase = :v)
    return pv,vl,vv,y
end

function bubble_temperature_impl(model,p,x,method::SingleFluidFlashMethod)
    Tl = _eval_inverse_generic_anc(model.ancillaries.fluid.saturation.data,p)
    vl = volume(model,p,Tl,phase = :l)
    vv = volume(model,p,Tl,phase = :v)
    return Tl,vl,vv,x
end

function dew_temperature_impl(model,p,y,method::SingleFluidFlashMethod)
    Tv = _eval_inverse_generic_anc(model.ancillaries.fluid.saturation.data2,p)
    vl = volume(model,p,Tv,phase = :l)
    vv = volume(model,p,Tv,phase = :v)
    return Tv,vl,vv,y
end

