"""
    SpecialComp(components::Vector{String},defaults=["water08"])

Auxiliary model that just stores the index of one component. faster than looking for the component string on each iteration. 

Used mainly for applying specific correlations in presence of a certain compound.
"""
struct SpecialComp <: EoSModel
    components::Vector{String}
    idx::Int
    references::Vector{String}
end

function SpecialComp(components::Vector{String},defaults=["water08"])
    idx = findfirst(in(defaults),components)
    idx === nothing && (idx = 0)
    return SpecialComp(components,idx,defaults)
end

function split_model(model::SpecialComp,subset=nothing)
    splitted = auto_split_model(model,subset)
    defaults = model.references
    for i in 1:length(splitted)
        splitted[i] = SpecialComp(splitted[i].components,defaults)
    end
    return splitted
end

Base.getindex(model::SpecialComp) = model.idx

