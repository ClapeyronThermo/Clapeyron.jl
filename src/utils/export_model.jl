# Export generic models
function export_model(model::EoSModel,name="")
    M = typeof(model)
    model_name = String(split(string(M),"{")[1])

    if hasfield(M,:groups)
        species = model.groups.flattenedgroups
        ncomps = length(species)
    else
        species = model.components
        ncomps = length(species)
    end


    like = Dict{Symbol,AbstractVector}()
    merge!(like,Dict(Symbol("species")=>species))

    species1 = Vector{String}()
    species2 = Vector{String}()

    for i in 1:ncomps-1
        append!(species1,fill(species[i],ncomps-i))
        append!(species2,species[i+1:end])
    end

    unlike = Dict{Symbol,AbstractVector}()
    merge!(unlike,Dict(Symbol("species1")=>species1))
    merge!(unlike,Dict(Symbol("species2")=>species2))

    assoc = Dict{Symbol,AbstractVector}()

    if hasfield(M,:params)
        P = typeof(model.params)
        params = fieldnames(P)
        for i in 1:length(params)
            if typeof(getfield(model.params,params[i])) <: SingleParam
                merge!(like,Dict(Symbol(params[i])=>getfield(model.params,params[i]).values))
            elseif typeof(getfield(model.params,params[i])) <: PairParam
                if params[i] == :sigma
                    merge!(like,Dict(Symbol(params[i])=>diagvalues(getfield(model.params,params[i]).values)*1e10))
                    binary = Vector{Float64}()
                    for j in 1:ncomps-1
                        append!(binary,getfield(model.params,params[i]).values[j+1:end,j]*1e10)
                    end
                    merge!(unlike,Dict(Symbol(params[i])=>binary))
                else
                    merge!(like,Dict(Symbol(params[i])=>diagvalues(getfield(model.params,params[i]).values)))
                    binary = Vector{Float64}()
                    for j in 1:ncomps-1
                        append!(binary,getfield(model.params,params[i]).values[j+1:end,j])
                    end
                    merge!(unlike,Dict(Symbol(params[i])=>binary))
                end
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

                merge!(assoc,Dict(Symbol(params[i])=>vals))
                if !any(keys(assoc).==:species1)
                    merge!(assoc,Dict(Symbol("species1")=>spe1))
                    merge!(assoc,Dict(Symbol("species2")=>spe2))
                    merge!(assoc,Dict(Symbol("site1")=>site1))
                    merge!(assoc,Dict(Symbol("site2")=>site2))
                end
            end
        end
    end

    if hasfield(M,:sites)
        site_types = model.sites.flattenedsites
        n_flatsites = model.sites.n_flattenedsites
        for i in 1:length(site_types)
            nsites = [n_flatsites[j][i] for j in 1:ncomps]
            merge!(like,Dict(Symbol("n_"*site_types[i])=>nsites))
        end
    end

    f = fieldnames(M)
    for i in f
        if typeof(getfield(model,i))<:EoSModel && hasfield(typeof(getfield(model,i)),:components) && i!=:vrmodel
            export_model(getfield(model,i),name)
        end
    end

    if name==""
        name = model_name
    else
        name=name*"_"*model_name
    end

    if length(like) != 1
        ParamTable(:like, Tables.columntable(like); name=name, location=".")
    end
    if length(unlike) != 2
        ParamTable(:unlike, Tables.columntable(unlike); name=name, location=".")
    end
    if !isempty(assoc)
        ParamTable(:assoc, Tables.columntable(assoc); name=name, location=".")
    end
end


# Export Activity coefficient models
function export_model(model::ActivityModel,name="")
    M = typeof(model)
    model_name = String(split(string(M),"{")[1])

    species = model.components
    ncomps = length(model)


    like = Dict{Symbol,AbstractVector}()
    merge!(like,Dict(Symbol("species")=>species))

    species1 = Vector{String}()
    species2 = Vector{String}()

    for i in 1:ncomps
        for j in 1:ncomps
            if i != j
                append!(species1,[species[i]])
                append!(species2,[species[j]])
            end
        end
    end

    unlike = Dict{Symbol,AbstractVector}()
    merge!(unlike,Dict(Symbol("species1")=>species1))
    merge!(unlike,Dict(Symbol("species2")=>species2))

    assoc = Dict{Symbol,AbstractVector}()

    if hasfield(M,:params)
        P = typeof(model.params)
        params = fieldnames(P)
        for i in 1:length(params)
            if typeof(getfield(model.params,params[i])) <: SingleParam
                merge!(like,Dict(Symbol(params[i])=>getfield(model.params,params[i]).values))
            elseif typeof(getfield(model.params,params[i])) <: PairParam
                binary = Vector{Float64}()
                for j in 1:ncomps
                    append!(binary,getfield(model.params,params[i]).values[1:end .!=j,j])
                end
                merge!(unlike,Dict(Symbol(params[i])=>binary))
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

                merge!(assoc,Dict(Symbol(params[i])=>vals))
                if !any(keys(assoc).==:species1)
                    merge!(assoc,Dict(Symbol("species1")=>spe1))
                    merge!(assoc,Dict(Symbol("species2")=>spe2))
                    merge!(assoc,Dict(Symbol("site1")=>site1))
                    merge!(assoc,Dict(Symbol("site2")=>site2))
                end
            end
        end
    end

    if hasfield(M,:sites)
        site_types = model.sites.flattenedsites
        n_flatsites = model.sites.n_flattenedsites
        for i in 1:length(site_types)
            nsites = [n_flatsites[j][i] for j in 1:ncomps]
            merge!(like,Dict(Symbol("n_"*site_types[i])=>nsites))
        end
    end

    f = fieldnames(M)
    for i in f
        if typeof(getfield(model,i))<:EoSModel && hasfield(typeof(getfield(model,i)),:components)
            export_model(getfield(model,i),name)
        end
    end

    if name==""
        name = model_name
    else
        name=name*"_"*model_name
    end

    if length(like) != 1
        ParamTable(:like, Tables.columntable(like); name=name, location=".")
    end
    if length(unlike) != 2
        ParamTable(:unlike, Tables.columntable(unlike); name=name, location=".")
    end
    if !isempty(assoc)
        ParamTable(:assoc, Tables.columntable(assoc); name=name, location=".")
    end
end
export export_model