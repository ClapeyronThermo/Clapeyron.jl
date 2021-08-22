abstract type MixingRule <:EoSModel end
function init_model(model::MixingRule,components,userlocations,verbose)
    return model
end

function init_model(model::Type{<:MixingRule},components,userlocations,verbose)
    verbose && @info("""Now creating alpha model:
    $idealmodel""")
    return model(components;userlocations,verbose)
end

has_sites(::Type{<:MixingRule})=false

include("vdW1f.jl")
include("Kay.jl")