# extract all sites
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
