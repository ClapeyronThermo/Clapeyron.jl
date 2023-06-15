
struct LinearMixing <: MixingRule end

is_splittable(::LinearMixing) = false

function LinearMixing(components;userlocations = String[],verbose = false)
    LinearMixing()
end

function v_scale(model::EmpiricMultiFluid,V,T,z,mixing::LinearMixing,∑z = sum(z))
    Vc = model.params.Tc.values
    dot(Vc,z)/∑z
end

function T_scale(model::EmpiricMultiFluid,V,T,z,mixing::LinearMixing,∑z = sum(z))
    Tc = model.params.Tc.values  
    return dot(Tc,z)/∑z
end
