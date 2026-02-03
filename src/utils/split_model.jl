#generated each_split_model for structs
function _each_split_model(param,field,fieldname,group,Ic,Ig)
    if !is_splittable(field)
        return field
    elseif isnothing(group)
        return each_split_model(field,Ic)
    else
        return each_split_model(field,group,Ic,Ig)
    end
end

function _each_split_model(param::EoSModel,field,fieldname,group,Ic,Ig)
    if !is_splittable(field) || fieldname == :references
        return field
    elseif fieldname == :components || isnothing(group) || field isa EoSModel
        return each_split_model(field,Ic)
    else
        return each_split_model(field,group,Ic,Ig)
    end
end

@generated function each_split_model_struct(param::P,group,Ic,Ig) where {P}
    all_fields = fieldnames(P)
    r = Expr(:call,Base.typename(P).wrapper)
    for field in all_fields
            field_sym = QuoteNode(field)
            push!(r.args,:(_each_split_model(param,param.$field,$field_sym,group,Ic,Ig)))
    end
    return r
end

each_split_model_struct(param,I) = each_split_model_struct(param,nothing,I,nothing)

each_split_model(model,group,I_component,I_group) = each_split_model(model,I_component)

function each_split_model(param::AbstractVector,I)
    val = param[I]
    eltype(param) <: AbstractArray && return deepcopy(val)
    return val
end

function each_split_model(param::UnitRange{Int},I)
    return 1:length(I)
end

each_split_model(param::ReferenceState,group,I_component,I_group) = each_split_model(param,I_component)

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

function each_split_model(param::PackedVofV,I)
    val = PackedVectorsOfVectors.pack(param[I])
    return val
end

function each_split_model(assoc::Compressed4DMatrix{T},I) where T
    len = length(assoc.values)
    iszero(len) && return Compressed4DMatrix{T}()
    old_idx = assoc.outer_indices
    idx_bool = findall(x -> (first(x) ∈ I) & (last(x) ∈ I),old_idx)
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

function each_split_model(param::ClapeyronParam,group,I_component,I_group)
    components = param.components
    if group === nothing
        return each_split_model(param,I_component)
    elseif components == group.components
        return each_split_model(param,I_component)
    elseif components == group.flattenedgroups
        return each_split_model(param,I_group)
    else
        __each_split_model_ambiguous_comps(param.name,typeof(param))
    end
end

@noinline function __each_split_model_ambiguous_comps(param,paramtype)
    throw(Argument("$param ($paramtype) is in a GC model, but does not have compatible component names for either component-based or group-based splitting."))
end

function each_split_model(param::EoSParam,group,Ic,Ig)
    each_split_model_struct(param,group,Ic,Ig)
end

function each_split_model(param::EoSParam,I)
    each_split_model_struct(param,I)
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
    if param.sites === nothing
        sites_i = nothing
    else
        sites_i = param.sites[I]
    end
    return AssocParam(
            param.name,
            param.components[I],
            _value,
            sites_i,
            param.sourcecsvs,
            param.sources
            )
end

function create_group_splitter(param::GroupParam,I)
    flattenedgroups = param.flattenedgroups
    len_groups = length(flattenedgroups)
    Ig = zeros(Int,len_groups)
    for i in I
        group_i = param.groups[i]
        for k in 1:length(group_i)
            j = findfirst(==(group_i[k]),flattenedgroups)::Int
            Ig[j] = j
        end
    end
    filter!(!iszero,Ig)
    return Ig
end

function each_split_model(param::GroupParam{TT},__group,Ic,Ig) where TT
    grouptype = param.grouptype
    components = param.components[Ic]
    groups = param.groups[Ic]
    n_groups = param.n_groups[Ic]
    sourcecsvs = param.sourcecsvs
    len_groups = length(param.flattenedgroups)

    flattenedgroups = param.flattenedgroups[Ig]
    i_groups = [[findfirst(isequal(group), flattenedgroups)::Int for group ∈ componentgroups] for componentgroups ∈ groups]
    n_flattenedgroups = Vector{Vector{TT}}(undef,length(Ic))

    #handling for intergroups
    n_intergroups = Vector{Matrix{TT}}(undef,length(Ic))
    empty_intergroup = Matrix{TT}(undef,(0,0))
    for (k,i) in pairs(Ic)
        pii = param.n_flattenedgroups[i]
        n_flattenedgroups[k] = pii[Ig]
        pij = param.n_intergroups[i]
        if !isempty(pij)
            n_intergroups[k] = pij[Ig,Ig]
        else
            n_intergroups[k] = empty_intergroup
        end
    end

    return GroupParam{TT}(
        components,
        groups,
        grouptype,
        n_groups,
        n_intergroups,
        i_groups,
        flattenedgroups,
        n_flattenedgroups,
        sourcecsvs)
end

function each_split_model(param::MixedGCSegmentParam{T},group,Ic,Ig) where T
    if length(param.values.v) == 0
        return MixedGCSegmentParam(param.name,param.components[Ic],deepcopy(param.values))
    end

    src = param.values
    ng = length(group.flattenedgroups)
    
    #count unique groups
    ncount = zeros(T,ng)
    for k in Ig
        ncount[k] = 1
    end
    ngg = count(!iszero,ncount)
    ncc = length(Ic)
    
    #reuse vector
    resize!(ncount,ngg*ncc)
    p = zeros(Int64,length(Ic)+1)
    p .= 1:ngg:(ncc*ngg + 1)
    dest = PackedVofV(p,ncount)
    
    for (k,i) in pairs(Ic)
        pii = src[i]
        true_n = @view(pii[Ig])
        dest[k] .= true_n
    end
    return MixedGCSegmentParam{T}(param.name,param.components[Ic],dest)
end

function each_split_model(group::GroupParam,I)
    Ig = create_group_splitter(group,I)
    return each_split_model(group,group,I,Ig)
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

function each_split_model(param::SiteParam,group,Ic,Ig)
    
    components = param.components
    if group === nothing
        site = each_split_model(param,Ic)
    elseif components == group.components
        site = each_split_model(param,Ic)
    elseif components == group.flattenedgroups
        site = each_split_model(param,Ig)
    else
        __each_split_model_ambiguous_comps("sites",SiteParam)
    end

    if group != nothing && site.site_translator != nothing
        ng = length(group.flattenedgroups)
        recalculate_site_translator!(site,Ig,ng)
    end
    
    return site
end

__split_site_translator(::Nothing,I) = nothing
__split_site_translator(s::Vector{Vector{NTuple{2,Int}}},I) = s[I]

function recalculate_site_translator!(sites::SiteParam,idxi,ng,bool_to_int = Int[])
    site_translator_i = sites.site_translator::Vector{Vector{NTuple{2,Int}}}
    #the tuple is (ki,site_kia) where ki is the position of the group, and site_ki is the site number in GC based sites
    #DO NOT use the second number. if you really need it, store the original SiteParam instead.
    #the first number is used to reference the gc pair at an specific component, via get_group_idx

    resize!(bool_to_int,ng)
    bool_to_int .= 0
    wk = 0
    for w in 1:ng
        if w in idxi
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

#EoSModel each_split_model code.

function each_split_model(model::EoSModel,I)
    if !is_splittable(model)
        return model
    end

    if I isa AbstractVector{Bool}
        return each_split_model(model,findall(I))
    end
    if has_groups(model)
        Ic = I
        groups = model.groups
        Ig = create_group_splitter(groups,I)
        return each_split_model_struct(model,groups,Ic,Ig)
    else
        return each_split_model_struct(model,I)
    end
end

"""
    split_model(model::EoSModel)
    split_model(model::EoSModel, splitter)

Takes in a model for a multi-component system and returns a vector of models. The result depends on the splitter used.

A model can be splitted in a list of submodels, where each submodel is of the same type as the original model, but it has a different combination of components.

Group-Contribution models are also splitted in a component basis, `split_model` takes care of converting between Group and Component basis automatically.

A splitter is just a list of indices for each submodel, valid splitters are:

- A list of integers: `split_model(model,[1,5,2])` will return three pure models.
- An integer: `split_model(model,1)` will return a list with one model corresponding to the first component
- A list of lists: `split_model(model,[[1,2],[3])` will return a list with two models, the first one will contain two components, the second one will be a pure model.

The default splitter is `1:length(model)`, that will return a list with all pure models.

## Examples

```julia-repl
julia> model = MonomerIdeal(["methane","propane","butane"])
MonomerIdeal with 3 components:
 "methane"
 "propane"
 "butane"
Contains parameters: Mw, reference_state

julia> split_model(model)
3-element Vector{MonomerIdeal}:
 MonomerIdeal("methane")
 MonomerIdeal("propane")
 MonomerIdeal("butane")
 
julia> split_model(model,[[1,2],[3,1]])
2-element Vector{MonomerIdeal}:
 MonomerIdeal("methane", "propane")
 MonomerIdeal("butane", "methane")

julia> split_model(model,2)
1-element Vector{MonomerIdeal}:
 MonomerIdeal("propane")
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

function split_model(param)
    if is_splittable(param)
        splitter = default_splitter(param)
        return split_model(param,splitter)
    else
        throw(ArgumentError("$param is not splittable, try passing an explicit splitter argument (`split_model(value,splitter)`)"))
    end
end

#general method
split_model(param,splitter) = _split_model(param,splitter)

function split_model(param::AbstractArray,splitter)
    s = size(param)
    length(s) > 1 && (@assert reduce(isequal,s))
    return _split_model(param,splitter)
end

#inner method to dispatch on type of splitter

#general splitter type is an iterator with eltype <: AbstractVector{Int}
function _split_model(param,splitter)
    if is_splittable(param) || param isa EoSModel
        return map(Base.Fix1(each_split_model,param),splitter)
    else
        return [fill(param,length(i)) for i ∈ splitter]
    end
end

function _split_model(param,splitter::AbstractVector{Int})
    if is_splittable(param) || param isa EoSModel
        f(i) = each_split_model(param,i:i)
        return map(f,splitter)
    else
        return [fill(param,length(i)) for i ∈ splitter]
    end
end

function _split_model(param,bool_splitter::AbstractVector{Bool})
    int_splitter = findall(bool_splitter)
    return _split_model(param,int_splitter)
end

_split_model(param,splitter::Nothing) = split_model(param)
_split_model(param,i::Int) = [each_split_model(param,i:i)]

for T in (:Symbol,:Tuple,:AbstractString,:Number,:Missing,:Nothing)
    @eval is_splittable(param::$T) = false
end

function _n_splitter(n)
    r = Vector{UnitRange{Int64}}(undef,n)
    for i in 1:n
        r[i] = i:i
    end
    return r
end

default_splitter(param::ClapeyronParam) = _n_splitter(length(param.components))
default_splitter(param::EoSModel) = 1:length(param)
default_splitter(param::AbstractArray) = _n_splitter(size(param,1))

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

"""

    split_pure_model(model,splitter)

Similar to `split_model` but promises that the result only has pure models. 
Some EoS models store a list of pure models and this function allows accessing that list.

"""
function split_pure_model(model)
    if is_splittable(model)
        splitter = default_splitter(model)
        return split_pure_model(model,splitter)
    else
        throw(ArgumentError("$model is not splittable, try passing an explicit splitter argument (`split_pure_model(value,splitter)`)"))
    end
end

split_pure_model(model,splitter) = split_model(model,splitter)
export split_model, split_model_binaries
