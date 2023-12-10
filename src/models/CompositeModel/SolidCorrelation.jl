"""
    SolidCorrelation{P,M,L} <: EoSModel

wrapper struct to signal that a `CompositeModel` uses solid correlations for the phase volume, melting and sublimation
"""
struct SolidCorrelation{P,M,L}  <: EoSModel
    components::Vector{String}
    phase::P
    melting::M
    sublimation::L
end
