abstract type LiquidVolumeModel <: EoSModel end
include("NaNLiquid/NaNLiquid.jl")
include("Rackett/RackettLiquid.jl")
include("COSTALD/costald.jl")
include("DIPPR105Liquid/DIPPR105Liquid.jl")
include("PolExpLiquid/PolExpLiquid.jl")
