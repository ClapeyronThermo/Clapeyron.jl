using Clapeyron, CSV, Tables
import Clapeyron: diagvalues

model = PR(["hexane"])

M = typeof(model)

species = model.components
ncomps = length(model)


like = Dict{String,Any}()
merge!(like,Dict("species"=>species))

species1 = Vector{String}()
species2 = Vector{String}()

for i in 1:ncomps-1
    append!(species1,fill(species[i],ncomps-i))
    append!(species2,species[i+1:end])
end

unlike = Dict{String,Any}()
merge!(unlike,Dict("species1"=>species1))
merge!(unlike,Dict("species2"=>species2))

assoc = Dict{String,Any}()

if hasfield(M,:params)
    P = typeof(model.params)
    params = fieldnames(P)
    for i in 1:length(params)
        if typeof(getfield(model.params,params[i])) <: SingleParam
            merge!(like,Dict(string(params[i])=>getfield(model.params,params[i]).values))
        elseif typeof(getfield(model.params,params[i])) <: PairParam
            merge!(like,Dict(string(params[i])=>diagvalues(getfield(model.params,params[i]).values)))
            binary = Vector{Float64}()
            for j in 1:ncomps-1
                append!(binary,getfield(model.params,params[i]).values[j+1:end,j])
            end
            merge!(unlike,Dict(string(params[i])=>binary))
        elseif typeof(getfield(model.params,params[i])) <: AssocParam
            site_types = model.sites.flattenedsites
            mat = getfield(model.params,params[i]).values
            nassoc = length(mat.values)
            spe1 = Vector{String}()
            spe2 = Vector{String}()
            site1 = Vector{String}()
            site2 = Vector{String}()
            vals = zeros(nassoc)
            for j in 1:nassoc
                outer_idx = mat.outer_indices[j]
                inner_idx = mat.inner_indices[j]
                push!(spe1,species[outer_idx[1]])
                push!(spe2,species[outer_idx[2]])
                push!(site1,site_types[inner_idx[1]])
                push!(site2,site_types[inner_idx[2]])
                vals[j] = mat.values[j]
            end

            merge!(assoc,Dict(string(params[i])=>vals))
            if !any(keys(assoc).==:species1)
                merge!(assoc,Dict("species1"=>spe1))
                merge!(assoc,Dict("species2"=>spe2))
                merge!(assoc,Dict("site1"=>site1))
                merge!(assoc,Dict("site2"=>site2))
            end
        end
    end
end

if hasfield(M,:sites)
    site_types = model.sites.flattenedsites
    n_flatsites = model.sites.n_flattenedsites
    for i in 1:length(site_types)
        nsites = [n_flatsites[j][i] for j in 1:ncomps]
        merge!(like,Dict("n_"*site_types[i]=>nsites))
    end
end

function export_csv(name,dict)
    header = collect(keys(dict))
    nparam = length(header)
    nentry = length(dict[header[1]])

    data = [dict[header[i]][j] for i in 1:nparam for j in 1:nentry]
    
    data = reshape(data,(nentry,nparam))
    df = Tables.table(data;header=header)
    CSV.write(name,df)
end

export_csv("like.csv",like)
export_csv("unlike.csv",unlike)
if !isempty(assoc)
    export_csv("assoc.csv",assoc)
end