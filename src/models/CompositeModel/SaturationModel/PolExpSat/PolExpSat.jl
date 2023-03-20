struct PolExpSat <: SaturationModel
    Tc::Float64
    Pc::Float64
    n::Vector{Float64}
    v::Vector{Float64}
end
#TODO: add Tmin here

function crit_pure(model::PolExpSat)
    return (model.Tc,model.Pc,NaN)
end

function saturation_pressure_impl(model::PolExpSat,T,method::SaturationCorrelation)
    nan = zero(T)/zero(T)
    Tc = model.Tc
    Pc = model.Pc
    T>Tc && return zero(T)/zero(T)
    Tr = T/Tc
    θ = 1.0-Tr
    lnPsatPc = evalexppoly(θ,model.n,model.v)/Tr
    Psat = exp(lnPsatPc)*Pc
    return Psat,nan,nan
end

Base.length(::PolExpSat) = 1
is_splittable(::PolExpSat) = false

export PolExpSat
