macro comps()
    return :($(esc(:(model.icomponents))))
end

macro groups()
    return :($(esc(:(model.iflattenedgroups))))
end

macro groups(component)
    return :($(esc(:(model.igroups[$(component)]))))
end

macro sites(component)
    return :($(esc(:(model.isites[$(component)]))))
end

macro f(func, args...)
    args = [esc(arg) for arg ∈ args]
    return :($(func)($(esc(:model)),$(esc(:V)),$(esc(:T)),$(esc(:z)),$(args...)))
end

macro newmodel(name, parent, paramstype)
    esc(quote struct $name{T <: IdealModel} <: $parent
        modelname::String
        components::Array{String,1}
        allcomponentsites::Array{Array{String,1},1}
        allcomponentnsites::Array{Array{Int,1},1}
        params::$paramstype
        icomponents::UnitRange{Int}
        isites::Array{UnitRange{Int},1}
        idealmodel::T
        references::Array{String,1}
        function $name(params::$paramstype, sites::SiteParam, idealmodel::T; references::Array{String,1}=String[]) where {T <: IdealModel}
            arbitraryparam = getfield(params, first(fieldnames($paramstype)))
            modelname = arbitraryparam.modelname
            components = arbitraryparam.components
            allcomponentsites = sites.allcomponentsites
            allcomponentnsites = sites.allcomponentnsites
            icomponents = 1:length(components)
            isites = [1:length(componentsites) for componentsites ∈ allcomponentsites]
            return new{T}(modelname, components, allcomponentsites, allcomponentnsites, params, icomponents, isites, idealmodel, references)
        end
        function $name(params::$paramstype, idealmodel::IdealModel; references::Array{String,1}=String[]) where {T <: IdealModel}
            arbitraryparam = getfield(params, first(fieldnames($paramstype)))
            components = arbitraryparam.components
            sites = SiteParam(components, [String[] for _ ∈ 1:length(arbitraryparam.components)], [Int[] for _ ∈ 1:length(arbitraryparam.components)], arbitraryparam.modelname)
            return $name(params, sites, idealmodel; references=references)
        end
        function $name(params::$paramstype, sites::SiteParam; references::Array{String,1}=String[]) where {T <: IdealModel}
            arbitraryparam = getfield(params, first(fieldnames($paramstype)))
            components = arbitraryparam.components
            idealmodel = BasicIdeal(components)
            return $name(params, sites, idealmodel; references=references)
        end
        function $name(params::$paramstype; references::Array{String,1}=String[]) where {T <: IdealModel}
            arbitraryparam = getfield(params, first(fieldnames($paramstype)))
            components = arbitraryparam.components
            println(components)
            idealmodel = BasicIdeal(components)
            sites = SiteParam(components, [String[] for _ ∈ 1:length(arbitraryparam.components)], [Int[] for _ ∈ 1:length(arbitraryparam.components)], arbitraryparam.modelname)
            return $name(params, sites, idealmodel; references=references)
        end
    end
    end)
end

macro newmodelgc(name, parent, paramstype)
    esc(quote struct $name{T <: IdealModel} <: $parent
        modelname::String
        components::Array{String,1}
        allcomponentgroups::Array{Array{String,1},1}
        allcomponentngroups::Array{Array{Int,1},1}
        flattenedgroups::Array{String,1}
        allcomponentnflattenedgroups::Array{Array{Int,1},1}
        allgroupsites::Array{Array{String,1},1}
        allgroupnsites::Array{Array{Int,1},1}
        params::$paramstype
        icomponents::UnitRange{Int}
        igroups::Array{Array{Int,1},1}
        iflattenedgroups::UnitRange{Int}
        isites::Array{UnitRange{Int},1}
        idealmodel::T
        references::Array{String,1}
        function $name(params::$paramstype, groups::GCParam, sites::SiteParam, idealmodel::IdealModel; references::Array{String,1}=String[]) where {T <: IdealModel}
            modelname = groups.modelname
            components = groups.components
            allcomponentgroups = groups.allcomponentgroups
            allcomponentngroups = groups.allcomponentngroups
            flattenedgroups = groups.flattenedgroups
            allcomponentnflattenedgroups = groups.allcomponentnflattenedgroups
            allgroupsites = sites.allcomponentsites
            allgroupnsites = sites.allcomponentnsites
            icomponents = 1:length(components)
            igroups = [[findfirst(isequal(group), flattenedgroups) for group ∈ componentgroups] for componentgroups ∈ allcomponentgroups]
            iflattenedgroups = 1:length(flattenedgroups)
            isites = [1:length(groupsites) for groupsites ∈ allgroupsites]
            return new{T}(modelname, components, allcomponentgroups, allcomponentngroups, flattenedgroups, allcomponentnflattenedgroups, allgroupsites, allgroupnsites, params, icomponents, igroups, iflattenedgroups, isites, idealmodel, references)
        end
        function $name(params::$paramstype, groups::GCParam, idealmodel::IdealModel; references::Array{String,1}=String[]) where {T <: IdealModel}
            modelname = groups.modelname
            components = groups.components
            sites = SiteParam(components, [String[] for _ ∈ 1:length(components)], [Int[] for _ ∈ 1:length(components)], modelname)
            return $name(params, groups, sites, idealmodel; references=references)
        end
        function $name(params::$paramstype, groups::GCParam, sites::SiteParam; references::Array{String,1}=String[]) where {T <: IdealModel}
            modelname = groups.modelname
            components = arbitraryparam.components
            idealmodel = BasicIdeal(components)
            return $name(params, groups, sites, idealmodel; references=references)
        end
        function $name(params::$paramstype, groups::GCParam; references::Array{String,1}=String[]) where {T <: IdealModel}
            modelname = groups.modelname
            components = groups.components
            idealmodel = BasicIdeal(components)
            sites = SiteParam(components, [String[] for _ ∈ 1:length(components)], [Int[] for _ ∈ 1:length(components)], modelname)
            return $name(params, groups, sites, idealmodel; references=references)
        end
    end
    end)
end
