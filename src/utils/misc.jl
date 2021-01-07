export create_z
function create_z(model::EoS, z::AbstractArray)
    return NamedArray(z, model.components)
end

# extract all sites
export extractsites
function extractsites(n_sites)
    return union([([collect(keys(n_sites[x])) for x in keys(n_sites)]...)...])
end

const LIQUID_STR = ("liquid","LIQUID","L","l")
is_liquid(str::String) = str in LIQUID_STR

const VAPOUR_STR = ("vapor","VAPOR","VAPOUR","vapour","g","G","v","V")
is_vapour(str::String) = str in VAPOUR_STR

const SUPERCRITICAL_STR = ("sc","SC","supercritical","SUPERCRITICAL")
is_supercritical(str::String) = str in SUPERCRITICAL_STR

function ∑(iterator)
    # wrapper for sum function that returns 0. if iterator is empty
    if isempty(collect(iterator))
        return 0.
    end 
    return sum(iterator)
end

"""
    eos_name(eos::EoS)::String

returns the name of the equation of state.
"""
eos_name(eos::EoS)::String = string(nameof(typeof(eos)))

