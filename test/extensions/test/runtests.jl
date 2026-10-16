using MultiComponentFlash: MultiComponentFlash
using Symbolics
using EoSSuperancillaries
using PythonCall
using Glenn
using Aqua
using Carnot

include("../../utils.jl")

#include_distributed("base/database.jl",4)
#include_distributed("base/solvers.jl",4)
#include_distributed("base/differentials.jl",4)
#include_distributed("base/misc.jl",4)

#include_distributed("base/qa.jl",5)
include("../../base/qa.jl")

#include_distributed("models/saft_pc.jl",1)
#include_distributed("models/cubic.jl",3)
#include_distributed("models/saft_others.jl",3)
#include_distributed("models/others.jl",1)
#include_distributed("models/saft_vr.jl",2)
#include_distributed("models/electrolytes.jl",1)

#include_distributed("methods/eos.jl",4)

#include_distributed("methods/eos_pcsaft.jl",5)
include("../../methods/eos_pcsaft.jl")

#include_distributed("methods/api.jl",2)
include("../../methods/mass_fractions.jl")
#include_distributed("methods/api_flash.jl",3)
#include_distributed("methods/electrolytes.jl",1)
#include_distributed("methods/estimation.jl",2)

#include_distributed("test_issues.jl",1)

include("extensions.jl")
include("carnot.jl")