abstract type TranslationModel <:EoSModel end
function init_model(model::TranslationModel,components,userlocations,verbose)
    return model
end

function init_model(model::Type{<:TranslationModel},components,userlocations,verbose)
    verbose && @info("""Now creating mixing model:
    $model""")
    return model(components;userlocations,verbose)
end

has_sites(::Type{<:TranslationModel})=false

include("NoTranslation.jl")
include("Rackett.jl")
include("Peneloux.jl")
include("MT.jl")