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
    file = CSV.File(DB_PATH*"/Activity/COSMOSAC/cosmo_inchikey.csv")
    CAS = file["CAS#"]
    INCHIKEY = file["INCHIKEY"]
    return CAS, INCHIKEY
end

cosmo_tol(_new,_old) = mapreduce((x,y) -> abs(x/y -1.0),+,_new,_old)

