include("consts.jl")
include("params.jl")
include("equations.jl")
include("GERG2008/GERG2008.jl")
include("EOS_LNG/EOS_LNG.jl")


struct EmpiricMultiFluidProperties <: EoSParam
    Mw::SingleParam{Float64}
end

struct EmpiricDepartureParam <: EoSParam end


struct EmpiricMultiFluid{𝔸} <: EmpiricHelmholtzModel
    components::Vector{String}
    params::EmpiricMultiFluidProperties
    pures::Vector{EmpiricSingleFluid{𝔸}}
    mixing::EmpiricDepartureParam
    references::Vector{String}
end

function EmpiricMultiFluid(components; pure_userlocations = String[],departure_userlocations = String[],verbose = false)
    pures = [SingleFluid(comp;userlocations = pure_userlocations,verbose = verbose) for comp in comps]
    Mw = SingleParam("Mw",components,[pure.properties.Mw for pure in pures])
    R0 = [pure.ideal.R0 for pure in pures]
    Rgas = [pure.properties.Rgas for pure in pures]
    references = reduce(vcat,pure.references for pure in pures)
end
