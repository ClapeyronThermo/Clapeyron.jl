struct ReactiveAqParams{T} <: EoSParam where T <: Float64
    ΔGf::SingleParam{T}
    ΔHf::SingleParam{T}
end

struct ReactiveAqModel{T}
    components::Vector{String}
    reactions::Vector
    params::ReactiveParams{Float64}
    eosmodel::T
end

function ReactiveAqModel(components::Vector{String}, reactions::Vector;
                       model::EoSModel,
                       userlocations=String[],
                       verbose=false)
    params = getparams(components, ["properties/formation_aqueous.csv"]; userlocations = userlocations, verbose = verbose)
    ΔGf = params["Gf"]
    ΔHf = params["Hf"]

    packagedparams = ReactiveParams(ΔGf,ΔHf)

    return ReactiveAqModel(components,reactions,packagedparams,model)
end

function stoichiometric_coefficient(model::ReactiveAqModel)
    ν = zeros(length(model.components),length(model.reactions))
    for (i,reaction) in enumerate(model.reactions)
        species = [reaction[k][1] for k in 1:length(reaction)]
        coeff = [reaction[k][2] for k in 1:length(reaction)]
        for (j,component) in enumerate(model.components)
            if component in species
                ν[j,i] = coeff[findfirst(isequal(component),species)]
            end
        end
    end
    return ν
end


function Base.show(io::IO, mime::MIME"text/plain", model::ReactiveAqModel)   
    print(io, typeof(model))
    length(model.eosmodel) == 1 && println(io, " with 1 component:")
    length(model.eosmodel) > 1 && println(io, " with ", length(model.eosmodel), " components:")
    show_pairs(io,model.components)
    println(io, "")

    length(model.reactions) == 1 && println(io, "and 1 reaction:")
    length(model.reactions) > 1 && println(io, "and ", length(model.reactions), " reactions:")
    for reaction in model.reactions
        species = [reaction[k][1] for k in 1:length(reaction)]
        coeff = [reaction[k][2] for k in 1:length(reaction)]
        reactants = ""
        products = ""
        for i in 1:length(species)
            if coeff[i] < 0
                reactants *= string(abs(coeff[i]))*" "*species[i]*" + "
            else
                products *= string(abs(coeff[i]))*" "*species[i]*" + "
            end
        end
        reactants = reactants[1:end-3]
        products = products[1:end-3]
        println(io, reactants, " <=> ", products)
    end
end

export ReactiveAqModel