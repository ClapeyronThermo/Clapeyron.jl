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
deferred_valtype(x) = Bool
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

@generated function recursive_deferred_valtype(x, ::Val{RECURSE}) where {RECURSE}
    # Build a nested promote_type(...) call entirely at specialisation time.
    T_acc = :(Bool)   # neutral starting point
    deferred_valtype_expr = Expr(:call,:(Base.promote_type))
    names = fieldnames(x)
    if length(names) == 0
        return :(Bool)
    end
    for (fname, ftype) in zip(fieldnames(x), fieldtypes(x))
        field = :(getfield(x, $(QuoteNode(fname))))

        ft = if ftype <: Number || ftype <: AbstractArray
            # Scalar numbers and arrays → grab element type
            :(eltype($field))
        elseif RECURSE
            # Recurse into nested structs
            :(recursive_deferred_valtype($field, Val{$RECURSE}()))
        else
            # Flat mode: delegate to the user-supplied hook
            :(deferred_valtype($field))
        end
        push!(deferred_valtype_expr.args,ft)
    end

    return deferred_valtype_expr   # returns a *type*, not a value
end
