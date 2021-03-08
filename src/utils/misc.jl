# extract all sites
const LIQUID_STR = (:liquid,:LIQUID,:L,:l)
is_liquid(sym::Symbol) = sym in LIQUID_STR
is_liquid(str::String) = is_liquid(Symbol(str)) 

const VAPOUR_STR = (:vapor,:VAPOR,:VAPOUR,:vapour,:g,:G,:v,:V,:gas,:GAS)
is_vapour(sym::Symbol) = sym in VAPOUR_STR
is_vapour(str::String) = is_vapour(Symbol(str))

const SUPERCRITICAL_STR = (:sc,:SC,:supercritical,:SUPERCRITICAL)
is_supercritical(sym::Symbol) = sym in SUPERCRITICAL_STR
is_supercritical(str::String) = is_vapour(Symbol(str))

function ∑(iterator) #not collecting is faster
    return reduce(Base.add_sum,iterator,init=0.0)
end