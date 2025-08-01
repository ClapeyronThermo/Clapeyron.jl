"""
    SolidCorrelation{P,M,L} <: RestrictedEquilibriaModel

Wrapper struct to signal that a `CompositeModel` uses solid correlations for the phase volume, melting and sublimation.
"""
struct SolidCorrelation{P,M,L} <: RestrictedEquilibriaModel
    components::Vector{String}
    phase::P
    melting::M
    sublimation::L
end

struct MeltingCorrelation <: ThermodynamicMethod end
