function auto_split_model(model::gcPCSAFTModel,subset=nothing)
        
    try
        allfields = Dict{Symbol,Any}()
        #if has_groups(typeof(model))
        #    raw_splitter = model.groups.i_groups
        #    subset !== nothing && throw("using subsets is not supported with Group Contribution models")
        #else
        #    raw_splitter = split_model(Vector(1:length(model.components)))
        #end

        comp_splitter = split_model(Vector(1:length(model.components)))

        if hasfield(typeof(model),:groups)
            gc_split = split_model(model.groups,comp_splitter)
            allfields[:groups] = gc_split
            allfields[:components] = split_model(model.groups.components,comp_splitter)
            gc_splitter = group_splitter(model.groups,gc_split)
        end

        len = length(comp_splitter)
        M = typeof(model)
        allfieldnames = fieldnames(M)

        #add here any special keys, that behave as non_splittable values
        for modelkey in [:references]
            if modelkey in allfieldnames
                if !haskey(allfields,modelkey)
                    allfields[modelkey] = fill(getproperty(model,modelkey),len)
                end
            end
        end

        for modelkey ∈ allfieldnames
            if !haskey(allfields,modelkey)
                modelx = getproperty(model,modelkey)
                if is_splittable(modelx)
                    if modelx isa EoSModel
                        allfields[modelkey]= split_model(modelx,subset)
                    elseif modelx isa SiteParam || modelx isa SiteTranslator
                        allfields[modelkey]= split_model(modelx,comp_splitter)
                    elseif modelx isa gcPCSAFTParam
                        allfields[modelkey]= auto_split_model(modelx,comp_splitter,gc_splitter)
                    else
                        allfields[modelkey]= split_model(modelx,gc_splitter)
                    end
                else
                    allfields[modelkey] = fill(modelx,len)
                end
            end
        end

        return [M((allfields[k][i] for k ∈ fieldnames(M))...) for i ∈ 1:len]::Vector{M}
    catch e
        M = typeof(model)
        @error "$M cannot be splitted"
        rethrow(e)
    end
end


function auto_split_model(model::gcPCSAFTParam,comp_splitter,gc_splitter)
    allfields = Dict{Symbol,Any}()

    len = length(comp_splitter)
    M = typeof(model)
    allfieldnames = fieldnames(M)

    for modelkey ∈ allfieldnames
        if !haskey(allfields,modelkey)
            modelx = getproperty(model,modelkey)
            if is_splittable(modelx)
                if modelx isa AssocParam
                    allfields[modelkey]= split_model(modelx,comp_splitter)
                else
                    allfields[modelkey]= split_model(modelx,gc_splitter)
                end
            else
                allfields[modelkey] = fill(modelx,len)
            end
        end
    end

    return [M((allfields[k][i] for k ∈ fieldnames(M))...) for i ∈ 1:len]::Vector{M}
end