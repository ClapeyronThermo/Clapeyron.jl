function each_split_model(param::AbstractVector,I)
    val = param[I]
    eltype(param) <: AbstractArray && return deepcopy(val)
    return val
end

function each_split_model(param::UnitRange{Int},I)
    return 1:length(I)
end

function each_split_model(param::ReferenceState,I)
    sym = param.std_type
    if length(param.a1) == 0
        return deepcopy(param)
    else
        return ReferenceState(param.components[I],param.a0[I],param.a1[I],param.T0,param.P0,param.H0[I],param.S0[I],param.z0[I],param.phase,param.std_type)
    end
end

function each_split_model(param::AbstractMatrix,I)
    val = param[I,I]
    eltype(param) <: AbstractArray && return deepcopy(val)
    return val
end

function each_split_model(y::SparseMatrixCSC{<:AbstractVector},I)
    x = y[I,I]
    m,n,colptr,rowval,nzval = x.m,x.n,x.colptr,x.rowval,x.nzval
    return SparseMatrixCSC(m,n,colptr,rowval,nzval)
end

function each_split_model(y::SparsePackedMofV,I)
    idx = y.idx[I,I]
    if iszero(length(y.storage))
        return SparsePackedMofV(y.storage,idx)
    end

    if iszero(nnz(idx))
        st = y.storage
        storage = PackedVofV([1],zeros(eltype(st.v),0))
        return SparsePackedMofV(storage,idx)
    else
        str = y.storage[nnz(idx)]
        storage = PackedVectorsOfVectors.pack(str)
        return SparsePackedMofV(storage,idx)
    end
end

function each_split_model(param::PackedVofV,I)
    val = PackedVectorsOfVectors.pack(param[I])
    return val
end

function each_split_model(assoc::Compressed4DMatrix{T},I) where T
    len = length(assoc.values)
    iszero(len) && return Compressed4DMatrix{T}()
    old_idx = assoc.outer_indices
    idx_bool = findall(x -> (first(x) ∈ I)&(last(x) ∈ I),old_idx)
    iszero(length(idx_bool)) && return Compressed4DMatrix{T}()
    values = assoc.values[idx_bool]
    outer_indices = assoc.outer_indices[idx_bool]
    inner_indices = assoc.inner_indices[idx_bool]
    out_val = length(I)
    outer_size = (out_val,out_val)
    inner_size = assoc.inner_size
    len2 = length(outer_indices)
    for i ∈ 1:len2
        i1,j1 = outer_indices[i]
        i2,j2 = findfirst(==(i1),I)::Int,findfirst(==(j1),I)::Int
        outer_indices[i] = (i2,j2)
    end
    return Compressed4DMatrix(values,outer_indices,inner_indices,outer_size,inner_size)
end

function each_split_model(param::SingleParameter,I)
    return SingleParameter(
        param.name,
        param.components[I],
        each_split_model(param.values,I),
        param.ismissingvalues[I],
        param.sourcecsvs,
        param.sources
    )
end

function each_split_model(param::PairParameter,I)
    value = each_split_model(param.values,I)
    ismissingvalues = param.ismissingvalues[I,I]
    components = param.components[I]
    res = PairParameter(
            param.name,
            components,
            value,
            ismissingvalues,
            param.sourcecsvs,
            param.sources
            )
    return res
end

function each_split_model(param::AssocParam,I)
    _value  = each_split_model(param.values,I)
    return AssocParam(
            param.name,
            param.components[I],
            _value,
            param.sites[I],
            param.sourcecsvs,
            param.sources
            )
end

function gc_each_split_model(param::GroupParam,I)
    grouptype = param.grouptype
    components = param.components[I]
    groups = param.groups[I]
    n_groups = param.n_groups[I]
    sourcecsvs = param.sourcecsvs

    #unique, but without allocating sets.
    _idx = zeros(Bool,length(param.flattenedgroups))
    for i in I
        group_i = param.groups[i]
        for k in 1:length(group_i)
            j = findfirst(==(group_i[k]),param.flattenedgroups)::Int
            _idx[j] = true
        end
    end

    len_groups = length(_idx)

    flattenedgroups = param.flattenedgroups[_idx]
    i_groups = [[findfirst(isequal(group), flattenedgroups)::Int for group ∈ componentgroups] for componentgroups ∈ groups]
    n_flattenedgroups = Vector{Vector{Int64}}(undef,length(I))
    for (k,i) in pairs(I)
        pii = param.n_flattenedgroups[i]
        n_flattenedgroups[k] = pii[_idx]
    end
    n_groups_cache  = PackedVectorsOfVectors.packed_fill(0.0,(length(ni) for ni in n_flattenedgroups))

    for (k,i) in pairs(I)
        pii = param.n_groups_cache[i]
        true_n = (pii[_idx])
        n_groups_cache[k] .= true_n
    end

    return _idx,GroupParam(
        components,
        groups,
        grouptype,
        n_groups,
        i_groups,
        flattenedgroups,
        n_flattenedgroups,
        n_groups_cache,
        sourcecsvs)
end

function gc_each_split_model(param::StructGroupParam,I)
    grouptype = param.grouptype
    components = param.components[I]
    groups = param.groups[I]
    n_groups = param.n_groups[I]
    sourcecsvs = param.sourcecsvs

    #unique, but without allocating sets.
    _idx = zeros(Bool,length(param.flattenedgroups))
    for i in I
        group_i = param.groups[i]
        for k in 1:length(group_i)
            j::Int = findfirst(==(group_i[k]),param.flattenedgroups)::Int
            _idx[j] = true
        end
    end

    len_groups = length(_idx)

    flattenedgroups = param.flattenedgroups[_idx]
    n_intergroups = Vector{Matrix{Int64}}(undef,length(I))
    i_groups = [[findfirst(isequal(group), flattenedgroups) for group ∈ componentgroups] for componentgroups ∈ groups]
    n_flattenedgroups = Vector{Vector{Int64}}(undef,length(I))
    for (k,i) in pairs(I)
        pii = param.n_flattenedgroups[i]
        n_flattenedgroups[k] = pii[_idx]

        pii = param.n_intergroups[i]
        n_intergroups[k] = pii[_idx,_idx]
    end
    n_groups_cache = PackedVectorsOfVectors.packed_fill(0.0,(length(ni) for ni in n_flattenedgroups))

    for (k,i) in pairs(I)
        pii = param.n_groups_cache[i]
        true_n = @view(pii[_idx])
        n_groups_cache[k] .= true_n
    end

    return _idx,StructGroupParam(
        components,
        groups,
        grouptype,
        n_groups,
        n_intergroups,
        i_groups,
        flattenedgroups,
        n_flattenedgroups,
        n_groups_cache,
        sourcecsvs)
end

function each_split_model(group::GroupParameter,I)
    _,gi = gc_each_split_model(group,I)
    return gi
end

function each_split_model(param::SiteParam,I)
    return SiteParam(
        param.components[I],
        param.sites[I],
        each_split_model(param.n_sites,I),
        param.i_sites[I],
        param.flattenedsites,
        param.n_flattenedsites[I],
        param.i_flattenedsites[I],
        param.sourcecsvs,
        __split_site_translator(param.site_translator,I))
end

__split_site_translator(::Nothing,I) = nothing
__split_site_translator(s::Vector{Vector{NTuple{2,Int}}},I) = s[I]


"""
    split_model(model::EoSModel)
Takes in a model for a multi-component system and returns a vector of model for each pure system.
## Example:
```julia-repl
julia> gerg2 = GERG2008(["propane","pentane"])
GERG008 model with 2 components:
"propane"
"pentane"
julia> split_model(gerg2)
2-element Vector{GERG2008}:
 GERG2008("propane")
 GERG2008("pentane")
```
"""
function split_model end

"""
    is_splittable(model)::Bool
Trait to determine if a `EoSModel` should be splitted by itself or can be simply filled into a vector.
This is useful in the case of models without any parameters, as those models are impossible by definition to split, because they don't have any underlying data.
The Default is `is_splittable(model) = true`.
"""
is_splittable(model) = true
is_splittable(null::Union{Nothing,Missing}) = false
is_splittable(::Number) = false
is_splittable(::AbstractString) = false
is_splittable(::Symbol) = false
is_splittable(::Tuple) = false


function split_model(param)
    if is_splittable(param)
        splitter = default_splitter(param)
        return split_model(param,splitter)
    else
        throw(ArgumentError("$param is not splittable, try passing an explicit splitter argument (`split_model(value,splitter)`)"))
    end
end


function split_model(param::ClapeyronParam,splitter)
    if is_splittable(param)
        return [each_split_model(param,i) for i ∈ splitter]
    else
        return [fill(param,length(i)) for i ∈ splitter]
    end
end

function split_model(param::AbstractArray,splitter)
    s = size(param)
    length(s) > 1 && (@assert reduce(isequal,s))
    return [each_split_model(param,i) for i ∈ splitter]
end

for T in (:Symbol,:Tuple,:AbstractString,:Number,:Missing,:Nothing)
    @eval split_model(param::$T,splitter) = [fill(param,length(i)) for i ∈ splitter]
end

function _n_splitter(n)
    r = Vector{UnitRange{Int64}}(undef,n)
    for i in 1:n
        r[i] = i:i
    end
    return r
end

default_splitter(param::ClapeyronParam) = _n_splitter(length(param.components))
default_splitter(param::EoSModel) = _n_splitter(length(param.components))
default_splitter(param::AbstractArray) = _n_splitter(size(param,1))

function split_model(Base.@nospecialize(params::EoSParam),splitter)
    T = typeof(params)
    split_paramsvals = (split_model(getfield(params,i),splitter) for i  ∈ fieldnames(T))
    return T.(split_paramsvals...)
end

#=
Specializations for splitting with groups
=#

function gc_eosparam_split_model(Base.@nospecialize(params::EoSParam),groups::GroupParameter,comp_splitter,gc_splitter)
    T = typeof(params)
    function _split(parami::ClapeyronParam)
        if parami.components == groups.components
            return split_model(parami,comp_splitter)
        elseif parami.components == groups.flattenedgroups
            return split_model(parami,gc_splitter)
        else
            throw(error("$parami is in a GC model, but does not have compatible component names for either component-based or group-based splitting."))
        end
    end
    _split(parami) = split_model(parami,gc_splitter)

    split_paramsvals = (_split(getfield(params,i)) for i  ∈ fieldnames(T))
    return T.(split_paramsvals...)
end

function group_splitter(group::T,splitter) where T <: GroupParameter
    n = length(splitter)
    group_split = Vector{T}(undef,n)
    idx_split = Vector{Vector{Bool}}(undef,n)
    flattenedgroups = group.flattenedgroups
    gc_splitter = Vector{Vector{Int}}(undef,n)
    m = length(group.flattenedgroups)
    for i in 1:n
        idxi,gi = gc_each_split_model(group,splitter[i])
        idx_split[i] = idxi
        group_split[i] = gi
        gc_splitter[i] = findall(isone,idxi)::Vector{Int}
    end
    return group_split,idx_split,gc_splitter
end

function recalculate_site_translator!(sites::Vector{SiteParam},idx_splitter)
    bool_to_int = Int[]
    for i in 1:length(idx_splitter)
        site_translator_i = sites[i].site_translator::Vector{Vector{NTuple{2,Int}}}
        #the tuple is (ki,site_kia) where ki is the position of the group, and site_ki is the site number in GC based sites
        #DO NOT use the second number. if you really need it, store the original SiteParam instead.
        #the first number is used to reference the gc pair at an specific component, via get_group_idx

        idxi = idx_splitter[i]
        resize!(bool_to_int,length(idxi))
        bool_to_int .= 0
        wk = 0
        for w in 1:length(idxi)
            if idxi[w]
                wk += 1
                bool_to_int[w] = wk
            end
        end
        for (l,s0) in pairs(site_translator_i)
            si = copy(s0)
            iszero(length(si)) && continue
            for a in 1:length(si)
                ki,_ = si[a]
                ki_new = bool_to_int[ki] #splitted group new index, if not zero
                si[a] = (ki_new,0)

            end
            site_translator_i[l] = filter!(x -> !iszero(first(x)),si)
        end
    end
end
#=
Start of EoSModel split_model functions
=#
split_model(model::EoSModel,splitter) = auto_split_model(model,splitter)

function auto_split_model(Base.@nospecialize(model::EoSModel),subset)
    try
        allfields = Dict{Symbol,Any}()

        M = typeof(model)
        allfieldnames = fieldnames(M)

        if subset === nothing
            splitter = _n_splitter(length(model.components))
        elseif eltype(subset) <: Integer
            splitter = [Int(i):Int(i) for i in subset]
        elseif eltype(subset) <: AbstractVector
            splitter = subset
        else
            throw(ArgumentError("Invalid type of subset. Expected subset::AbstractVector{Union{Int,AbstractVector{Int}}} or subset::Nothing"))
        end

        len = length(splitter)
        _has_groups = has_groups(M)
        #shortcut for directly non-splittable models

        if _has_groups
            gc_split,idx_splitter,gc_splitter = group_splitter(model.groups,splitter)
            allfields[:groups] = gc_split
            allfields[:components] = split_model(model.groups.components::Vector{String},splitter)
            comp_splitter = splitter
            splitter = gc_splitter
        else
            comp_splitter = splitter
            gc_splitter = splitter
        end
        #add here any special keys, that behave as non_splittable values
        for modelkey in (:references,)
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
                    if modelx isa SiteParam && _has_groups && modelx.site_translator !== nothing
                        #process site_translator
                        split_sites = split_model(modelx,comp_splitter)
                        recalculate_site_translator!(split_sites,idx_splitter)
                        allfields[modelkey] = split_sites
                    elseif modelx isa ClapeyronParam && _has_groups
                        #in this particular case, we can suppose that we have the components field
                        if modelx.components == model.groups.flattenedgroups
                            allfields[modelkey] = split_model(modelx,gc_splitter)
                        elseif modelx.components == model.groups.components
                            allfields[modelkey] = split_model(modelx,comp_splitter)
                        else
                            throw(error("$modelx is in a GC model, but does not have compatible component names for either component-based or group-based splitting."))
                        end
                    elseif modelx isa EoSParam && _has_groups
                        #we supppose a EoSParam has only one layer of splitting
                        allfields[modelkey] = gc_eosparam_split_model(modelx,model.groups,comp_splitter,gc_splitter)
                    else
                        #we suppose that this can splitted on his own (EoSModels are here.)
                        allfields[modelkey] = split_model(modelx,comp_splitter)
                    end
                else
                    allfields[modelkey] = fill(modelx,len)
                end
            end
        end

        return [M((allfields[k][i] for k ∈ allfieldnames)...) for i ∈ 1:len]::Vector{M}
    catch e
        M = typeof(model)
        @error "$M cannot be splitted"
        rethrow(e)
    end
end

##fallback,around 50 times slower if there is any need to read csvs

function simple_split_model(Base.@nospecialize(model::EoSModel),subset = nothing)
    MODEL = typeof(model)
    pure = Vector{MODEL}(undef,0)
    if subset === nothing
        comps = model.components
    else
        comps = model.components[subset]
    end
    for comp ∈ comps
        push!(pure,MODEL([comp]))
    end
    return pure
end

"""
    split_model_binaries(model::EoSModel)::Vector{EoSModel}

Given a multicomponent `EoSModel`, returns a list with the combination of all binary models.

## Example
```julia-repl
julia> model = PCPSAFT(["water", "methanol", "propyleneglycol","methyloxirane"])
PCPSAFT{BasicIdeal} with 4 components:
 "water"
 "methanol"
 "propyleneglycol"
 "methyloxirane"
Contains parameters: Mw, segment, sigma, epsilon, dipole, dipole2, epsilon_assoc, bondvol

julia> split_model_binaries(model)
6-element Vector{PCPSAFT{BasicIdeal}}:
 PCPSAFT{BasicIdeal}("water", "methanol")
 PCPSAFT{BasicIdeal}("water", "propyleneglycol")
 PCPSAFT{BasicIdeal}("water", "methyloxirane")
 PCPSAFT{BasicIdeal}("methanol", "propyleneglycol")
 PCPSAFT{BasicIdeal}("methanol", "methyloxirane")
 PCPSAFT{BasicIdeal}("propyleneglycol", "methyloxirane")
```
"""
function split_model_binaries(model)
    idx = Vector{Int}[]
    n = length(model)
    for i in 1:n
        for j in i+1:n
            push!(idx,[i,j])
        end
    end
    split_model(model,idx)
end

export split_model, split_model_binaries
