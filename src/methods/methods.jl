molecular_weight(model::EoSModel,z = @SVector [1.]) = 0.001*mapreduce(+,*,paramvals(model.params.Mw),z)
include("initial_guess.jl")
include("differentials.jl")
include("vt.jl")
include("property_solvers.jl")
include("pt.jl")
#include("unitful.jl")


