struct PolExpSat <: SaturationModel
    Tc::Float64
    Pc::Float64
    n::Vector{Float64}
    v::Vector{Float64}
end

function crit_pure(model::PolExpSat)
    return (model.Tc,model.Pc,NaN)
end

function saturation_pressure_impl(model::PolExpSat,T,method::SaturationCorrelation)
    nan = zero(T)/zero(T)
    Tc = model.Pc
    Pc = model.Tc
    T>T_c && return zero(T)/zero(T)
    Tr = T/Tc
    θ = 1.0-Tr
    lnPsatPc = evalexppoly(θ,model.n,model.v)
    Psat = exp(lnPsatPc)*Pc
    return Psat
    return psat,nan,nan
end

Base.length(::PolExpSat) = false
is_splittable(::PolExpSat) = false
