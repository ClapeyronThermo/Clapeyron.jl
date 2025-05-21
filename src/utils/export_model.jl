using OrderedCollections

Base.@nospecialize

"""
    export_model(model::EoSModel,name="";location=".")
Exports model parameters to CSVs. Unless the `name` kwarg is specified, the name of the files will follow the convention `singledata_EoS`, `pairdata_EoS` and `assocdata_EoS`. Files will be saved within the current directory unless the `location` argument is specified.

Note that it will export all submodel parameters (e.g. Alpha function parameters for cubic EoS).
"""
function export_model(model::EoSModel,name="";location=".")
    M = typeof(model)
    model_name = summary(model)

    if hasfield(M,:groups)
        species = model.groups.flattenedgroups
        ncomps = length(species)
    else
        species = model.components
        ncomps = length(species)
    end

    if hasfield(M,:params)
        P = typeof(model.params)
        params = fieldnames(P)
        export_like(model,params,name,location,species,ncomps)
        export_unlike(model,params,name,location,species,ncomps)
        export_assoc(model,params,name,location,species,ncomps)
    end


    f = fieldnames(M)
    for i in f
        fi = getfield(model,i)
        if getfield(model,i) isa EoSModel && hasfield(fieldtype(M,i),:components) && i != :vrmodel
            export_model(fi,name;location=location)
        end
    end
end

export export_model

function export_like(model::EoSModel,params,name,location,species,ncomps)
    M = typeof(model)
    P = typeof(model.params)
    model_name = summary(model)

    like = OrderedDict{Symbol,AbstractVector}()
    like[:species] = species

    for i in 1:length(params)
        paramtype = fieldtype(P,i)
        paramname = params[i]
        paramvalue = getfield(model.params,params[i])
        if paramtype <: SingleParam
            like[paramname] = paramvalue.values
        elseif paramtype <: PairParam
            if paramname == :sigma
                like[:sigma] = diagvalues(paramvalue.values) * 1e10
            elseif all(!iszero,diagvalues(paramvalue.values)) #all nonzero diagonal values
                like[paramname] = diagvalues(paramvalue.values)
            end
        end
    end

    if has_sites(model)
        sites = getsites(model)
        site_types = sites.flattenedsites
        n_flatsites = sites.n_flattenedsites
        for i in 1:length(site_types)
            nsites = [n_flatsites[j][i] for j in 1:ncomps]
            like[Symbol("n_"*site_types[i])] = nsites
        end
    end

    if name==""
        name = model_name
    else
        name=name*"_"*model_name
    end

    if length(like) != 1
        ParamTable(:like, Tables.columntable(like); name=name, location=location)
    end
end

function export_unlike(model::EoSModel,params,name,location,species,ncomps)
    M = typeof(model)
    P = typeof(params)
    model_name = summary(model)

    species1 = Vector{String}()
    species2 = Vector{String}()

    for i in 1:ncomps-1
        append!(species1,fill(species[i],ncomps-i))
        append!(species2,species[i+1:end])
    end

    unlike = OrderedDict{Symbol,AbstractVector}()
    unlike[:species1] = species1
    unlike[:species2] = species2

    for i in 1:length(params)
        paramtype = fieldtype(P,i)
        paramname = params[i]
        paramvalue = getfield(model.params,params[i])

        if paramvalue isa PairParameter
            binary = Vector{Float64}()
            if params[i] == :sigma
                for j in 1:ncomps-1
                    append!(binary,paramvalue.values[j+1:end,j]*1e10)
                end
            else
                for j in 1:ncomps-1
                    append!(binary,paramvalue.values[j+1:end,j])
                end
            end
            unlike[paramname] = binary
        end
    end

    if name==""
        name = model_name
    else
        name=name*"_"*model_name
    end

    if length(unlike) != 2
        ParamTable(:unlike, Tables.columntable(unlike); name=name, location=location)
    end
end

function export_unlike(model::ActivityModel,params,name,location,species,ncomps)
    M = typeof(model)
    P = typeof(model.params)
    model_name = summary(model)

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

    unlike = OrderedDict{Symbol,AbstractVector}()
    unlike[:species1] = species1
    unlike[:species2] = species2

    for i in 1:length(params)
        paramname = params[i]
        paramvalue = getfield(model.params,params[i])
        if paramvalue isa PairParameter
            binary = Vector{Float64}()
            for j in 1:ncomps
                append!(binary,paramvalue.values[j,1:end .!=j])
            end
            unlike[paramname] = binary
        end
    end

    if name==""
        name = model_name
    else
        name=name*"_"*model_name
    end

    if length(unlike) != 2
        ParamTable(:unlike, Tables.columntable(unlike); name=name, location=location)
    end
end

function export_unlike(model::ABCubicModel,params,name,location,species,ncomps)
    M = typeof(model)
    P = typeof(model.params)
    model_name = summary(model)

    species1 = Vector{String}()
    species2 = Vector{String}()

    for i in 1:ncomps-1
        append!(species1,fill(species[i],ncomps-i))
        append!(species2,species[i+1:end])
    end

    unlike = OrderedDict{Symbol,AbstractVector}()
    unlike[:species1] = species1
    unlike[:species2] = species2

    for i in 1:length(params)
        paramname = params[i]
        paramvalue = getfield(model.params,params[i])
        if paramvalue isa PairParameter
            if params[i] == :a
                binary = Vector{Float64}()
                for j in 1:ncomps-1
                    aij = paramvalue.values[j+1:end,j]
                    aj = paramvalue.values[j,j]
                    ai = diagvalues(paramvalue.values)[j+1:end]
                    kij = @. 1-aij/(sqrt(ai*aj))
                    append!(binary,kij)
                end
                unlike[:k] = binary
            elseif params[i] == :b
                binary = Vector{Float64}()
                for j in 1:ncomps-1
                    bij = paramvalue.values[j+1:end,j]
                    bj = paramvalue.values[j,j]
                    bi = diagvalues(paramvalue.values)[j+1:end]
                    lij = @. 1-2*bij/(bi+bj)
                    append!(binary,lij)
                end
                unlike[:l] = binary
            end
        end
    end

    if name==""
        name = model_name
    else
        name=name*"_"*model_name
    end
    @show unlike

    if length(unlike) != 2
        ParamTable(:unlike, Tables.columntable(unlike); name=name, location=location)
    end
end

function export_assoc(model::EoSModel,params,name,location,species,ncomps)
    M = typeof(model)
    P = typeof(model.params)
    model_name = summary(model)

    assoc = OrderedDict{Symbol,AbstractVector}()

    for i in 1:length(params)
        paramname = params[i]
        paramvalue = getfield(model.params,params[i])
        if paramvalue isa AssocParam
            site_types = model.sites.flattenedsites
            mat = paramvalue.values
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

            if !any(keys(assoc).==:species1)
                assoc[:species1] = spe1
                assoc[:species2] = spe2
                assoc[:site1] = site1
                assoc[:site2] = site2
            end

            assoc[paramname] = vals
        end
    end

    if name==""
        name = model_name
    else
        name=name*"_"*model_name
    end

    if !isempty(assoc)
        ParamTable(:assoc, Tables.columntable(assoc); name=name, location=location)
    end
end

Base.@specialize