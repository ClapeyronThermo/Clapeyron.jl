abstract type EoSModel end

abstract type SAFTModel <: EoSModel end
abstract type GCSAFTModel <: SAFTModel end
abstract type CubicModel <: EoSModel end
abstract type IdealModel <: EoSModel end
abstract type EmpiricHelmholtzModel <: EoSModel end
abstract type EoSParam end

function getsites(pairs::Dict{String,SingleParam{Int}})
    arbitraryparam = first(values(pairs))
    components = arbitraryparam.components
    allcomponentsites = arbitraryparam.allcomponentsites
    sourcecsvs = unique([([x.sourcecsvs for x in values(pairs)]...)...])
    allcomponentnsites = [[pairs[allcomponentsites[i][j]].values[i] for j ∈ 1:length(allcomponentsites[i])] for i ∈ 1:length(components)]  # or groupsites
    return SiteParam(components, allcomponentsites, allcomponentnsites, sourcecsvs)
end

function idealmodelselector(idealmodelstring::String, components::Array{String,1}; verbose::Bool=false)
    normalisedidealmodelstring = normalisestring(idealmodelstring)
    if normalisedidealmodelstring == "" 
        modelsym = :basic
    else
        modelsym = Symbol(normalisedidealmodelstring)
    end

    return idealmodelselector(Val{modelsym},components,verbose=verbose)
end

function idealmodelselector(::Type{Val{:monomer}},components;verbose=false)
    return MonomerIdeal(components; verbose=verbose)
end

function idealmodelselector(::Type{Val{:reid}},components;verbose=false)
    return ReidIdeal(components; verbose=verbose)
end

function idealmodelselector(::Type{Val{:walker}},components;verbose=false)
    return WalkerIdeal(components; verbose=verbose)
end

function idealmodelselector(::Type{Val{:basic}},components;verbose=false)
    return BasicIdeal(components; verbose=verbose)
end

function getsites(pairs::Dict{String,SingleParam{Int}})
    arbitraryparam = first(values(pairs))
    components = arbitraryparam.components
    allcomponentsites = arbitraryparam.allcomponentsites
    sourcecsvs = unique([([x.sourcecsvs for x in values(pairs)]...)...])
    allcomponentnsites = [[pairs[allcomponentsites[i][j]].values[i] for j ∈ 1:length(allcomponentsites[i])] for i ∈ 1:length(components)]  # or groupsites
    return SiteParam(components, allcomponentsites, allcomponentnsites, sourcecsvs)
end

function idealmodelselector(idealmodelstring::String, components::Array{String,1}; verbose::Bool=false)
    #=normalisedidealmodelstring = normalisestring(idealmodelstring)
    if normalisedidealmodelstring == "monomer"
        return MonomerIdeal(components; verbose=verbose)
    elseif normalisedidealmodelstring == "reid"
        return ReidIdeal(components; verbose=verbose)
    elseif normalisedidealmodelstring == "walker"
        return WalkerIdeal(components; verbose=verbose)
    elseif normalisedidealmodelstring == "basic" || normalisedidealmodelstring == ""
        return BasicIdeal(components; verbose=verbose)
    else
        error("Your selected ideal model ", idealmodelstring, " is not recognised.")
    end
    =#
    normalisedidealmodelstring = normalisestring(idealmodelstring)
    if normalisedidealmodelstring == "" 
        modelsym = :basic
    else
        modelsym = Symbol(normalisedidealmodelstring)
    end

    return idealmodelselector(Val{modelsym},components,verbose=verbose)
end

function idealmodelselector(::Type{Val{:monomer}},components;verbose=false)
    return MonomerIdeal(components; verbose=verbose)
end

function idealmodelselector(idealmodel::Type{T}, components::Array{String,1}; verbose::Bool=false) where T <: IdealModel
    return idealmodel(components,verbose=verbose)
end