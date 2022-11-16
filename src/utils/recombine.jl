function recombine!(model::EoSModel)
    if has_sites(model)
        recombine!(model.sites)
    end
    
    if has_groups(model)
        recombine!(model.groups)
    end

    ideal = idealmodel(model)
    
    if ideal !== nothing
        recombine!(ideal)
    end

    recombine_impl!(model) #this has to be user defined.
    return model
end

