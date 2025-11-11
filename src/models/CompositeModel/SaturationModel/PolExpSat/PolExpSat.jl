struct PolExpSat <: SaturationModel
    data::GenericAncEvaluator
    data2::GenericAncEvaluator #used in pseudo pure models.
end

PolExpSat(data::GenericAncEvaluator) = PolExpSat(data,data)

PolExpSat(anc::Solvers.ChebyshevRange) = PolExpSat(GenericAncEvaluator(anc))

function crit_pure(model::PolExpSat)
    single_component_check(crit_pure,model)
    if model.data.type == :superanc
        Tc = model.data.superanc.range[end]
        Pc = _eval_generic_anc(model.data,Tc)
        return Tc, Pc, NaN
    end
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
