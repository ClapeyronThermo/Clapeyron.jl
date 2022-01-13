function init_model(model::AlphaModel,components,userlocations,verbose)
    return model
end

function init_model(model::Type{<:AlphaModel},components,userlocations,verbose)
    verbose && @info("""Now creating alpha model:
    $model""")
    return model(components;userlocations,verbose)
end

include("NoAlpha.jl")
include("RKAlpha.jl")
include("PRAlpha.jl")
include("CPAAlpha.jl")
include("sCPAAlpha.jl")
include("PR78Alpha.jl")
include("soave.jl")
include("BM.jl")
include("Twu.jl")
include("MT.jl")