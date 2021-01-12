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
        lengthcomponents::Int
        icomponents::UnitRange{Int}

        allcomponentsites::Array{Array{String,1},1}
        lengthallcomponentsites::Array{Int,1}
        allcomponentnsites::Array{Array{Int,1},1}
        isites::Array{UnitRange{Int},1}

        params::$paramstype
        idealmodel::T
        absolutetolerance::Float64
        references::Array{String,1}

        function $name(params::$paramstype, sites::SiteParam, idealmodel::T;
                       references::Array{String,1}=String[],
                       absolutetolerance::Float64=1E-12) where {T <: IdealModel}
            arbitraryparam = getfield(params, first(fieldnames($paramstype)))
            modelname = arbitraryparam.modelname

            components = arbitraryparam.components
            lengthcomponents = length(components)
            icomponents = 1:lengthcomponents

            allcomponentsites = sites.allcomponentsites
            lengthallcomponentsites = [length(componentsites) for componentsites ∈ allcomponentsites]
            allcomponentnsites = sites.allcomponentnsites
            isites = [1:lengthallcomponentsites[i] for i ∈ icomponents]

            return new{T}(modelname, components, lengthcomponents, icomponents,
                          allcomponentsites, lengthallcomponentsites, allcomponentnsites, isites,
                          params, idealmodel, absolutetolerance, references)
        end
        function $name(params::$paramstype, idealmodel::T;
                       references::Array{String,1}=String[],
                       absolutetolerance::Float64=1E-12) where {T <: IdealModel}
            arbitraryparam = getfield(params, first(fieldnames($paramstype)))
            components = arbitraryparam.components
            sites = SiteParam(components, [String[] for _ ∈ 1:length(arbitraryparam.components)],
                [Int[] for _ ∈ 1:length(arbitraryparam.components)], arbitraryparam.modelname)
            return $name(params, sites, idealmodel; references=references)
        end
        function $name(params::$paramstype, sites::SiteParam;
                       references::Array{String,1}=String[],
                       absolutetolerance::Float64=1E-12)
            arbitraryparam = getfield(params, first(fieldnames($paramstype)))
            components = arbitraryparam.components
            idealmodel = BasicIdeal(components)
            return $name(params, sites, idealmodel; references=references)
        end
        function $name(params::$paramstype;
                       references::Array{String,1}=String[],
                       absolutetolerance::Float64=1E-12)
            arbitraryparam = getfield(params, first(fieldnames($paramstype)))
            components = arbitraryparam.components
            idealmodel = BasicIdeal(components)
            sites = SiteParam(components, [String[] for _ ∈ 1:length(arbitraryparam.components)],
                [Int[] for _ ∈ 1:length(arbitraryparam.components)], arbitraryparam.modelname)
            return $name(params, sites, idealmodel; references=references)
        end
    end
    end)
end

macro newmodelgc(name, parent, paramstype)
    esc(quote struct $name{T <: IdealModel} <: $parent
        modelname::String

        components::Array{String,1}
        lengthcomponents::Int
        icomponents::UnitRange{Int}

        allcomponentgroups::Array{Array{String,1},1}
        lengthallcomponentgroups::Array{Int,1}
        allcomponentngroups::Array{Array{Int,1},1}
        igroups::Array{Array{Int,1},1}

        flattenedgroups::Array{String,1}
        lengthflattenedgroups::Int
        allcomponentnflattenedgroups::Array{Array{Int,1},1}
        iflattenedgroups::UnitRange{Int}

        allgroupsites::Array{Array{String,1},1}
        lengthallgroupsites::Array{Int,1}
        allgroupnsites::Array{Array{Int,1},1}
        isites::Array{UnitRange{Int},1}

        params::$paramstype
        idealmodel::T
        absolutetolerance::Float64
        references::Array{String,1}

        function $name(params::$paramstype, groups::GCParam, sites::SiteParam, idealmodel::T;
                       references::Array{String,1}=String[],
                       absolutetolerance::Float64=1E-12) where {T <: IdealModel}
            modelname = groups.modelname

            components = groups.components
            lengthcomponents = length(components)
            icomponents = 1:lengthcomponents

            flattenedgroups = groups.flattenedgroups
            lengthflattenedgroups = length(flattenedgroups)
            allcomponentnflattenedgroups = groups.allcomponentnflattenedgroups
            iflattenedgroups = 1:lengthflattenedgroups

            allcomponentgroups = groups.allcomponentgroups
            lengthallcomponentgroups = [length(allcomponentgroups[i]) for i in icomponents]
            allcomponentngroups = groups.allcomponentngroups
            igroups = [[findfirst(isequal(group), flattenedgroups) for group ∈ componentgroups] for componentgroups ∈ allcomponentgroups]

            allgroupsites = sites.allcomponentsites
            lengthallgroupsites = [length(groupsites) for groupsites ∈ allgroupsites]
            allgroupnsites = sites.allcomponentnsites
            isites = [1:lengthallgroupsites[k] for k ∈ iflattenedgroups]

            return new{T}(modelname, components, lengthcomponents, icomponents,
                          allcomponentgroups, lengthallcomponentgroups, allcomponentngroups, igroups,
                          flattenedgroups, lengthflattenedgroups, allcomponentnflattenedgroups, iflattenedgroups,
                          allgroupsites, lengthallgroupsites, allgroupnsites, isites,
                          params, idealmodel, absolutetolerance, references)
        end
        function $name(params::$paramstype, groups::GCParam, idealmodel::T;
                       references::Array{String,1}=String[],
                       absolutetolerance::Float64=1E-12) where {T <: IdealModel}
            modelname = groups.modelname
            components = groups.components
            sites = SiteParam(components, [String[] for _ ∈ 1:length(components)], [Int[] for _ ∈ 1:length(components)], modelname)
            return $name(params, groups, sites, idealmodel; references=references, absolutetolerance=absolutetolerance)
        end
        function $name(params::$paramstype, groups::GCParam, sites::SiteParam;
                       references::Array{String,1}=String[],
                       absolutetolerance::Float64=1E-12) 
            modelname = groups.modelname
            components = groups.components
            idealmodel = BasicIdeal(components)
            return $name(params, groups, sites, idealmodel; references=references, absolutetolerance=absolutetolerance)
        end
        function $name(params::$paramstype, groups::GCParam;
                       references::Array{String,1}=String[],
                       absolutetolerance::Float64=1E-12) 
            modelname = groups.modelname
            components = groups.components
            idealmodel = BasicIdeal(components)
            sites = SiteParam(components, [String[] for _ ∈ 1:length(components)], [Int[] for _ ∈ 1:length(components)], modelname)
            return $name(params, groups, sites, idealmodel; references=references, absolutetolerance=absolutetolerance)
        end
    end
    end)
end

