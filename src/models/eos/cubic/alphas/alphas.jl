function init_model(model::AlphaModel,components,userlocations,verbose)
    return model
end

function init_model(model::Type{<:AlphaModel},components,userlocations,verbose)
    verbose && @info("""Now creating alpha model:
    $model""")
    return model(components;userlocations,verbose)
end

has_sites(::Type{<:AlphaModel})=false
has_groups(::Type{<:AlphaModel})=false

include("RKAlpha.jl")
include("PRAlpha.jl")
include("PR78Alpha.jl")
include("soave.jl")
include("BM.jl")
include("Twu.jl")
include("MT.jl")