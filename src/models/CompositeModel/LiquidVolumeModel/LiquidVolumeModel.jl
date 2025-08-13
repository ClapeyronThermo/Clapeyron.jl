abstract type LiquidVolumeModel <: GibbsBasedModel end



include("NaNLiquid/NaNLiquid.jl")
include("Rackett/RackettLiquid.jl")
include("COSTALD/costald.jl")
include("DIPPR105Liquid/DIPPR105Liquid.jl")
include("PolExpLiquid/PolExpLiquid.jl")
include("GrenkeElliott/GrenkeElliottWater.jl")
include("HoltenWater/HoltenWater.jl")