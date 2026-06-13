#preferred_valtype

if isdefined(Base,:Fix)
    @eval begin
        deferred_valtype(obj::Base.Fix{N,F,X}) where {N,F,X} = deferred_valtype(obj.x)
    end
else
    deferred_valtype(obj::Base.Fix1{F,X}) where {F,X} = deferred_valtype(obj.x)
    deferred_valtype(obj::Base.Fix2{F,X}) where {F,X} = deferred_valtype(obj.x)
end

deferred_valtype(x::T) where T <: Number = T
deferred_valtype(x::AbstractArray{T}) where T <: Number = T
deferred_valtype(::Type{T}) where T <: Number = T
deferred_valtype(::Type{AbstractArray{T}}) where T <: Number = T

#inner_function

inner_function(f) = f

if isdefined(Base,:Fix)
    @eval begin
        inner_function(obj::Base.Fix{N,F,X}) where {N,F,X} = inner_function(obj.f)
    end
else
    inner_function(obj::Base.Fix1{F,X}) where {F,X} = inner_function(obj.x)
    inner_function(obj::Base.Fix2{F,X}) where {F,X} = inner_function(obj.x)
end

function auto_context(f::F,tag::TT = ∂Tag{inner_function(f)}(),recursive::Val{RECURSIVE} = Val{false}()) where {F,TT,RECURSIVE}
    V = recursive_deferred_valtype(f,recursive)
    return WithContext{TT,V,F}(f)
end

skip_recursive_deferred_valtype(::Type{Nothing}) = true
skip_recursive_deferred_valtype(::Type{Missing}) = true


@generated function recursive_deferred_valtype(f::F,::Val{RECURSIVE}) where {F,RECURSIVE}
    names = fieldnames(F)
    types = fieldtypes(F)
    deferred_valtype_expr = Expr(:call,:(Base.promote_type))
    for (name,type) in zip(names,types)
        skip_recursive_deferred_valtype(type) && continue
        if RECURSIVE
            if hasmethod(deferred_valtype,Tuple{type})
                push!(deferred_valtype_expr.args,:(deferred_valtype(x.$name)))
            else
                push!(deferred_valtype_expr.args,:(recursive_deferred_valtype(x.$name,Val{true}())))
            end
        else
            push!(deferred_valtype_expr.args,:(deferred_valtype(x.$name)))
        end
    end
    return deferred_valtype_expr
end
