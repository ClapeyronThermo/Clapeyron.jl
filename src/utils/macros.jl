macro comps()
    return :($(esc(:(model.components))))
end

macro groups(component)
    return :($(esc(:(keys(model.group_multiplicities[$(component)])))))
end

macro groups()
    return :($(esc(:(model.groups))))
end

macro sites(component)
    return :($(esc(:(keys(model.params.n_sites[$(component)])))))
end

# macro sites()
#     return :($(esc(:(model.sites))))
# end

macro f(func, args...)
    args = [esc(arg) for arg ∈ args]
    return :($(func)($(esc(:model)),$(esc(:z)),$(esc(:V)),$(esc(:T)),$(args...)))
end

macro newmodel(name, parent, paramstype)
    esc(quote struct $name <: $parent
        modelname::String
        components::Array{String,1}
        allcomponentsites::Array{Array{String,1},1}
        allncomponentsites::Array{Array{Int,1},1}
        params::$paramstype
        icomponents::UnitRange{Int}
        isites::Array{UnitRange{Int},1}
        function $name(params::$paramstype, sites::SiteParam)
            arbitraryparam = getfield(params, first(fieldnames($paramstype)))
            modelname = arbitraryparam.modelname
            components = arbitraryparam.components
            allcomponentsites = sites.allcomponentsites
            allncomponentsites = sites.allncomponentsites
            icomponents = 1:length(components)
            isites = [1:length(componentsites) for componentsites ∈ allcomponentsites]
            return new(modelname, components, allcomponentsites, allncomponentsites, params, icomponents,isites)
        end
        function $name(params::$paramstype)
            arbitraryparam = getfield(params, first(fieldnames($paramstype)))
            sites = SiteParam(arbitraryparam.components, [String[] for _ ∈ 1:length(arbitraryparam.components)], [Int[] for _ ∈ 1:length(arbitraryparam.components)], arbitraryparam.modelname)
            return $name(params, sites)
        end
    end
    end)
end

macro newmodelgc(name, parent, paramstype)
    esc(quote struct $name <: $parent
        modelname::String
        components::Array{String,1}
        allcomponentgroups::Array{Array{String,1},1}
        allncomponentgroups::Array{Array{Int,1},1}
        flattenedgroups::Array{String,1}
        allcomponentsites::Array{Array{String,1},1}
        allncomponentsites::Array{Array{Int,1},1}
        params::$paramstype
        icomponents::UnitRange{Int}
        icomponentgroups::Array{Array{Int,1},1}
        isites::Array{UnitRange{Int},1}
        function $name(params::$paramstype, groups::GCParam, sites::SiteParam)
            modelname = groups.modelname
            components = groups.components
            allcomponentgroups = groups.allcomponentgroups
            allncomponentgroups = groups.allncomponentgroups
            flattenedgroups = groups.flattenedgroups
            allcomponentsites = sites.allcomponentsites
            allncomponentsites = sites.allncomponentsites
            icomponents = 1:length(components)
            icomponentgroups = [[findfirst(isequal(group), flattenedgroups) for group ∈ componentgroups] for componentgroups ∈ allcomponentgroups]
            isites = [1:length(componentsites) for componentsites ∈ allcomponentsites]
            return new(modelname, components, allcomponentgroups, allncomponentgroups, flattenedgroups, allcomponentsites, allncomponentsites, params, icomponents, icomponentgroups, isites)
        end
        function $name(params::$paramstype, groups::GCParam)
            sites = SiteParam(groups.components, [String[] for _ ∈ 1:length(groups.components)], [Int[] for _ ∈ 1:length(groups.components)], groups.modelname)
            return $name(params, groups, sites)
        end
    end
    end)
end
