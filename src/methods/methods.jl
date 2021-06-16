molecular_weight(model::EoSModel,z = @SVector [1.]) = 0.001*mapreduce(+,*,mw(model),z)
mw(model::EoSModel) = paramvals(model.params.Mw)

include("initial_guess.jl")
include("differentials.jl")
include("VT.jl")
include("property_solvers/property_solvers.jl")
include("pT.jl")
include("unitful_base.jl")
include("unitful_methods.jl")

