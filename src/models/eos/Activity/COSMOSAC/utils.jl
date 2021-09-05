function COSMO_parse_Pi(param::SingleParam{String})
    Vec = Vector{Vector{Float64}}(undef,0)
    n = length(param.components)
    default = Float64[]
    for i in 1:n
        if !param.ismissingvalues[i]
            push!(Vec,parse.(Float64,split(param.values[i])))
        else
            push!(Vec,default)
        end
    end
    SingleParam(param.name,param.components,Vec,param.ismissingvalues,param.sourcecsvs,param.sources)
end

cosmo_tol(_new,_old) = mapreduce((x,y) -> abs(x/y -1.0),+,_new,_old)

