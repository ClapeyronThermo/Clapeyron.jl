"""
    recombine_mixing!(model::MultiFluid,mixing::MixingRule)

Calculates missing values, using the parameters stored in `model`. modifies `mixing` implace. this function is called at `MultiFluid` model creation time.

"""
recombine_mixing!(model::MultiFluid,mixing,estimate) = nothing

function init_multifluid_mixing(model::EoSModel,components,estimate,userlocations = String[],verbose = false)
    return model
end

function init_multifluid_mixing(::Type{𝕄},components,estimate,userlocations = String[],verbose = false) where 𝕄
    if verbose
        if estimate == :off
            additional = "no estimation"
        elseif estimate == :lb
            additional = "Lorentz-Berthelot estimation"
        elseif estimate == :linear
            additional = "linear estimation" 
        end
        @info "Building an instance of $(info_color(string(𝕄))) with components $components, using $additional of missing mixing parameters"
    end
    return 𝕄(components;userlocations,estimate,verbose)
end

include("Asymmetric.jl")
include("LB.jl")
include("linear.jl")
include("TL.jl")

v_scale(model::MultiFluid,z,mixing::Nothing,∑z = sum(z)) = v_scale(model,z,LinearMixing(),∑z)
T_scale(model::MultiFluid,z,mixing::Nothing,∑z = sum(z)) = T_scale(model,z,LinearMixing(),∑z)

