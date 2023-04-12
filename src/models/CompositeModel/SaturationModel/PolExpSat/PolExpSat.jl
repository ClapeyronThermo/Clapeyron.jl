struct PolExpSat <: SaturationModel
    data::GenericAncEvaluator
end
  
function crit_pure(model::PolExpSat)
    return (model.data.input_r,model.output_r,NaN)
end

function saturation_pressure_impl(model::PolExpSat,T,method::SaturationCorrelation)
    nan = zero(T)/zero(T)

    return Psat,nan,nan
end

Base.length(::PolExpSat) = 1
is_splittable(::PolExpSat) = false

export PolExpSat
