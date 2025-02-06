function recombine!(model::EoSModel)
    #normally non-splittable models doesn't have component-dependent parameters, they are "final" in a way.
    if !is_splittable(model)
        return model
    end
    
    if has_sites(model)
        recombine!(getsites(model))
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
