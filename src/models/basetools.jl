abstract type EoS end

abstract type SAFTModel <: EoS end
abstract type CubicModel <: EoS end

abstract type EoSParam end

function getsites(pairs::Dict{String,SingleParam{Int}})
    arbitraryparam = first(values(pairs))
    components = arbitraryparam.components
    allcomponentsites = arbitraryparam.allcomponentsites
    modelname = arbitraryparam.modelname
    allncomponentsites = [[pairs[allcomponentsites[i][j]].values[i] for j ∈ 1:length(allcomponentsites[i])] for i ∈ 1:length(components)]
    return SiteParam(components, allcomponentsites, allncomponentsites, modelname)
end
