
#=
Liquid Cp models are functions of only temperature.

the main difference with a IdealModel is that there is no xlogx term and no log(V) term

=#
abstract type LiquidCpModel <: EoSModel end

function volume_impl(model::LiquidCpModel,p,T,z,phase,threaded,vol0)
    return zero(Base.promote_eltype(model,p,T,z))
end

function reference_state(model::LiquidCpModel)
    return __reference_state(model)
end

idealmodel(::LiquidCpModel) = ZeroIdeal() #this is a hack to exclude the ideal part.

include("PolynomialCpLiquid/PolynomialCpLiquid.jl")
include("PolynomialCpLiquid/ConstantCpLiquid.jl")




