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

function promote_model(::Type{T},model::EoSModel) where T <: Number
    return promote_model_struct(T,model)
end

function promote_model(::Type{T},model::Array) where T <: Number
    return promote_model.(T,model)
end

@generated function promote_model_struct(::Type{T},model::M) where {T,M}
    names = fieldnames(M)
    Base.typename(M).wrapper
    expr = Expr(:call,Base.typename(M).wrapper)
    M̄ = parameterless_type(M)
    for name in names
        if isconcretetype(fieldtype(M̄,name))
            push!(expr.args,:(deepcopy(model.$name))) #TODO: this should only copy in some particular cases.
        else
            push!(expr.args,:(promote_model($T,model.$name)))
        end
    end
    return expr
end

function promote_model(::Type{T},param::SingleParam) where T
    values = Vector{T}(param.values)
    return SingleParam{T}(param.name,param.components,values,param.ismissingvalues,param.sourcecsvs,param.sources)
end

function promote_model(::Type{T},param::PairParam) where T
    values = Matrix{T}(param.values)
    return PairParam{T}(param.name,param.components,values,param.ismissingvalues,param.sourcecsvs,param.sources)
end

function promote_model(::Type{T},param::AssocParam) where T
    v = param.values
    vals = Vector{T}(v.values)
    values = Compressed4DMatrix(vals,v.outer_indices,v.inner_indices,v.outer_size,v.inner_size)
    return AssocParam{T}(param.name,param.components,values,param.sites,param.sourcecsvs,param.sources)
end

function promote_model(::Type{T},param::MixedGCSegmentParam) where T
    vals = param.values
    p,v = vals.p,vals.v
    v2 = Vector{T}(v)
    values = PackedVofV(p,v2)
    return MixedGCSegmentParam{T}(param.name,param.components,values)
end
promote_model(::Type{T},v::AbstractArray{<:AbstractString}) where T = v

promote_model(::Type{T},params::EoSParam) where T = promote_model_struct(T,params)
promote_model(::Type{T},param::ClapeyronParam) where T = deepcopy(param)
promote_model(::Type{T},param::AssocOptions) where T = param
promote_model(::Type{T},param::ReferenceState) where T = param
promote_model(::Type{T},param::SpecialComp) where T = param
promote_model(::Type{T},param::GroupParam) where T = param
promote_model(::Type{T},param::SiteParam) where T = param

