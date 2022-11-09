function recombine!(model::EoSModel)
    if hassites(model)
        recombine!(model.sites)
    end
    
    if hasgroups(model)
        recombine!(model.groups)
    end

    recombine_impl!(model) #this has to be user defined.
    return model
end

