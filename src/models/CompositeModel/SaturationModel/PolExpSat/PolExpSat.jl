struct PolExpSat <: SaturationModel
    data::GenericAncEvaluator
end
  
function crit_pure(model::PolExpSat)
    single_component_check(crit_pure,model)
    anc = model.data
    return (anc.input_r,anc.output_r,NaN)
end

function saturation_pressure_impl(model::PolExpSat,T,method::SaturationCorrelation)
    nan = zero(T)/zero(T)
    Psat = _eval_generic_anc(model.data,T)
    return Psat,nan,nan
end

Base.length(::PolExpSat) = 1
is_splittable(::PolExpSat) = false

export PolExpSat
