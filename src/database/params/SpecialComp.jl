"""
    SpecialComp(components::Vector{String},defaults=["water08"])

Auxiliary model that just stores the index of one component. faster than looking for the component string on each iteration. 

Used mainly for applying specific correlations in presence of a certain compound.
"""
struct SpecialComp <: ClapeyronParam
    components::Vector{String}
    idx::Int
    defaults::Vector{String}
end

function SpecialComp(components::Vector{String},defaults=["water08"])
    idx = findfirst(in(defaults),components)
    idx === nothing && (idx = 0)
    return SpecialComp(components,idx,defaults)
end

function each_split_model(param::SpecialComp,I)
    components = param.components[I]
    defaults = param.defaults
    idx = findfirst(in(defaults),components)
    idx === nothing && (idx = 0)
    return SpecialComp(components,idx,copy(defaults))
end

Base.getindex(model::SpecialComp) = model.idx
#for compatibility with parametric eos param:
Base.eltype(model::SpecialComp) = Bool
