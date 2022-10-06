abstract type LiquidVolumeModel <: EoSModel end

include("Rackett/RackettLiquid.jl")
include("COSTALD/costald.jl")
include("DIPPR105Liquid/DIPPR105Liquid.jl")

