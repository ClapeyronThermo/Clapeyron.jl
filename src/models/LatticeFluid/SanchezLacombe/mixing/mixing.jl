abstract type SLMixingRule end

function init_model(model::SLMixingRule,components,userlocations,verbose)
    return model
end

function init_model(model::Type{<:SLMixingRule},components,userlocations,verbose)
    verbose && @info("""Now creating Sanchez-Lacombe mixing rule:
    $model""")
    return model(components;userlocations,verbose)
end
