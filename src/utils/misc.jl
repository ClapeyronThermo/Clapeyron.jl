# extract all sites
const LIQUID_STR = ("liquid","LIQUID","L","l")
is_liquid(str::String) = str in LIQUID_STR

const VAPOUR_STR = ("vapor","VAPOR","VAPOUR","vapour","g","G","v","V")
is_vapour(str::String) = str in VAPOUR_STR

const SUPERCRITICAL_STR = ("sc","SC","supercritical","SUPERCRITICAL")
is_supercritical(str::String) = str in SUPERCRITICAL_STR

function ∑(iterator) #not collecting is faster
    return reduce(Base.add_sum,iterator,init=0.0)
end
