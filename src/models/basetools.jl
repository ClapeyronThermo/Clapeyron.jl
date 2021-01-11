abstract type EoSModel end

abstract type SAFTModel <: EoSModel end
abstract type CubicModel <: EoSModel end
abstract type IdealModel <: EoSModel end

abstract type EoSParam end

function getsites(pairs::Dict{String,SingleParam{Int}})
    arbitraryparam = first(values(pairs))
    components = arbitraryparam.components
    allcomponentsites = arbitraryparam.allcomponentsites
    modelname = arbitraryparam.modelname
    allcomponentnsites = [[pairs[allcomponentsites[i][j]].values[i] for j ∈ 1:length(allcomponentsites[i])] for i ∈ 1:length(components)]  # or groupsites
    return SiteParam(components, allcomponentsites, allcomponentnsites, modelname)
end

function idealmodelselector(idealmodelstring::String, components::Array{String,1})
    normalisedidealmodelstring = normalisestring(idealmodelstring)
    if normalisedidealmodelstring == "monomer"
        return MonomerIdeal(components)
    elseif normalisedidealmodelstring == "reid"
        return ReidIdeal(components)
    elseif normalisedidealmodelstring == "walker"
        return WalkerIdeal(components)
    elseif normalisedidealmodelstring == "basic" || normalisedidealmodelstring == ""
        return BasicIdeal(components)
    else
        error("Your selected ideal model ", idealmodelstring, " is not recognised.")
    end
end
