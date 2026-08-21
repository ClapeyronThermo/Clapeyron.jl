using Test
t1 = @elapsed using Clapeyron
using CoolProp #CoolProp ext

ENV["JULIA_TEST_FAILFAST"] = "false"
ENV["CLAPEYRON_SHOW_REFERENCES"] = "FALSE" #Disable showing citations

@info "Loading Clapeyron took $(round(t1,digits = 2)) seconds"
@info "Coolprop: $(Clapeyron.is_coolprop_loaded())"

if true #by default, we autodetect if we are on CI or not, use `if false` to eliminate CI autodetection 
    include("ci.jl")
else
    __run_test(value) = true
end

function include_distributed(path,value::Int) #function to distribute workers in CI
    if __run_test(value)
        Base.include(@__MODULE__(), path)
    end
end

#file with test utilities
include("utils.jl")

#actual files
include_distributed("base/database.jl",4)
include_distributed("base/solvers.jl",4)
include_distributed("base/differentials.jl",4)
include_distributed("base/misc.jl",4)
include_distributed("base/qa.jl",5)

include_distributed("models/saft_pc.jl",1)
include_distributed("models/cubic.jl",3)
include_distributed("models/saft_others.jl",3)
include_distributed("models/others.jl",1)
include_distributed("models/saft_vr.jl",2)
include_distributed("models/electrolytes.jl",1)

include_distributed("methods/eos.jl",4)
include_distributed("methods/eos_pcsaft.jl",5)
include_distributed("methods/api.jl",2)
include_distributed("methods/api_flash.jl",3)
include_distributed("methods/electrolytes.jl",1)
include_distributed("methods/estimation.jl",2)

include_distributed("test_issues.jl",1)
