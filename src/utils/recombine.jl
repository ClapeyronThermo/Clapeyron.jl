function recombine!(model::EoSModel)
    if has_sites(model)
        recombine!(model.sites)
    end
    
    if has_groups(model)
        recombine!(model.groups)
    end

    if hasfield(typeof(model),:idealmodel)
        ideal = idealmodel(model)
        if ideal !== nothing
            recombine!(ideal)
        end
    end

    recombine_impl!(model) #this has to be user defined.
    return model
end

function recombine_impl!(model::EoSModel)
    return model
end
