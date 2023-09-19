function COSMO_parse_Pi(param::SingleParam{String})
    Vec = Vector{Vector{Float64}}(undef,0)
    n = length(param.components)
    
    
    default = Float64[]
    for i in 1:n
        if !param.ismissingvalues[i]
            push!(Vec,_vecparser(param.values[i]))
        else
            push!(Vec,default)
        end
    end
    SingleParam(param.name,param.components,Vec,param.ismissingvalues,param.sourcecsvs,param.sources)
end

function get_cosmo_comps()
    file = String(take!(Downloads.download("https://raw.githubusercontent.com/usnistgov/COSMOSAC/master/profiles/UD/complist.txt",IOBuffer())))
    lines = split(file,r"\n")
    words = split.(lines," ")
    CAS = [words[i][3] for i in 2:2263]
    INCHIKEY = [words[i][7] for i in 2:2263]
    return CAS, INCHIKEY
end

cosmo_tol(_new,_old) = mapreduce((x,y) -> abs(x/y -1.0),+,_new,_old)

