module ClapeyronMeasurementsExt
using Clapeyron
using Clapeyron.ForwardDiff
using Measurements

using Measurements: Measurement
using ForwardDiff: Dual

#if promote_rule(::Type{Measurement},Type{D}) is defined, then we suppose the overloads are already loaded.
include("ClapeyronMeasurementsExt/rules.jl")

end