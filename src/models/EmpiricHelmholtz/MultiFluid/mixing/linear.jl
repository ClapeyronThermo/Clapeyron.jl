
struct LinearMixing <: MixingRule end

is_splittable(::LinearMixing) = false

function LinearMixing(components;userlocations = String[],verbose = false)
    LinearMixing()
end

function v_scale(model::MultiFluid,z,mixing::LinearMixing,∑z)
    Vc = model.params.Tc.values
    dot(Vc,z)/∑z
end

function T_scale(model::MultiFluid,z,mixing::LinearMixing,∑z)
    Tc = model.params.Tc.values  
    return dot(Tc,z)/∑z
end

export LinearMixing