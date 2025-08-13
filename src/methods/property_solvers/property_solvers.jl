#volume calculations
include("volume.jl")

#lnϕ calculations
include("fugacity_coefficient.jl")

#single component properties and solvers
include("singlecomponent/singlecomponent.jl")

#multiple component properties and solvers
include("multicomponent/multicomponent.jl")

#stability calculations
include("stability/tpd.jl")
include("stability/stability.jl")

#T-X and P-X solvers
include("Tproperty.jl")
include("Pproperty.jl")

#spinodal solvers
include("spinodal.jl")

