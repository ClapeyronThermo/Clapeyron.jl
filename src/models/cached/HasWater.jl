const WATER_STRINGS = String["water"]

"""
    HasWater(components::Vector{String})

Auxiliary model that just stores the index of the water component. faster than looking for the water string on each iteration. 

If you need to add another component that is water, but it isn't called water, like `water-laffite2016`, just push that string to `Clapeyron.WATER_STRINGS`

"""
struct HasWater <: EoSModel
    components::Vector{String}
    idx::Int
end

function HasWater(components::Vector{String})
    idx = findfirst(in(WATER_STRINGS),components)
    idx === nothing && (idx = 0)
    return HasWater(components,idx)
end

function split_model(model::HasWater,subset=nothing)
    splitted = auto_split_model(model,subset)
    for i in 1:length(splitted)
        splitted[i] = HasWater(splitted[i].components)
    end
    return splitted
end

Base.getindex(model::HasWater) = model.idx

