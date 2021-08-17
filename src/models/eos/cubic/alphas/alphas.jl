function initialize_idealmodel(model::AlphaModel,components,userlocations,verbose)
    return model
end

function initialize_idealmodel(model::Type{<:AlphaModel},components,userlocations,verbose)
    verbose && @info("""Now creating alpha model:
    $idealmodel""")
    return model(components;userlocations,verbose)
end


include("soave.jl")