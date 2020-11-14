export create_z
function create_z(model::EoS, z::AbstractArray)
    return NamedArray(z, model.components)
end

export gc
function gc(args::Union{Tuple{String, Array{Pair{String, Int64},1}}, String, Tuple{String}}...)
    """
    The format for arguments is either the species name (if it available in the OpenSAFT database,
    or a tuple consisting of the species name, followed by a list of group => multiplicity pairs.
    For example:
    @GC( "sp1",
        ("sp2", ["CH3" => 2, "CH2" => 4]),
        ("sp3", ["CH3" => 2, "CH2" => 6]))
    """
    # This will be generalised soon
    species = Dict()
    for arg in args
        if typeof(arg) <: Tuple{String, Array{Pair{String, Int64},1}}
            string_key_dict = Dict(arg[2])
            species[Set([arg[1]])] = DefaultDict(0, Dict(Set([key]) => value for (key, value) in string_key_dict))
        elseif typeof(arg) <: Union{String, Tuple{String}}
            filepath = joinpath(dirname(pathof(OpenSAFT)), "../database", "SAFTgammaMie", "groups.csv")
            foundlines = findmatches(filepath, arg)
            if isempty(foundlines)
                error(arg * " not found.")
            end
            string = parseline(filepath, foundlines[end])[2]
            string_key_dict = Dict(eval(Meta.parse(string)))
            species[Set([arg])] = DefaultDict(0, Dict(Set([key]) => value for (key, value) in string_key_dict))
        end
    end
    return species
end

# extract all sites
export extractsites
function extractsites(n_sites)
    return union([([collect(keys(n_sites[x])) for x in keys(n_sites)]...)...])
end

function ∑(iterator)
    # wrapper for sum function that returns 0. if iterator is empty
    if isempty(collect(iterator))
        return 0.
    end 
    return sum(iterator)
end

