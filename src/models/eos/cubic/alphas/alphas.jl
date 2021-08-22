function init_model(model::AlphaModel,components,userlocations,verbose)
    return model
end

function init_model(model::Type{<:AlphaModel},components,userlocations,verbose)
    verbose && @info("""Now creating alpha model:
    $idealmodel""")
    return model(components;userlocations,verbose)
end

include("RKAlpha.jl")
include("PRAlpha.jl")
include("soave.jl")
include("BM.jl")