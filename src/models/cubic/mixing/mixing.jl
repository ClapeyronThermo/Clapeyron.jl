abstract type MixingRule <:EoSModel end
function init_model(model::MixingRule,components,activity,userlocations,activity_userlocations,verbose)
    return model
end

function init_model(model::Type{<:MixingRule},components,activity,userlocations,activity_userlocations,verbose)
    verbose && @info("""Now creating mixing model:
    $model""")
    return model(components;activity,userlocations,activity_userlocations,verbose)
end

include("vdW1f.jl")
include("Kay.jl")
include("HV.jl")
include("MHV1.jl")
include("MHV2.jl")
include("LCVM.jl")
include("WS.jl")
include("PSRK.jl")
include("VTPR.jl")
include("UMR.jl")