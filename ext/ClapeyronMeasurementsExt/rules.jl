using ForwardDiff: Dual, DiffRules, NaNMath, LogExpFunctions, SpecialFunctions
using Measurements: Measurement
import Base: +,-,/,*,promote_rule 

function promote_rule(::Type{Measurement{V}}, ::Type{Dual{T, V, N}}) where {T,V,N}
    Dual{Measurement{T}, V, N}
end

function promote_rule(::Type{Measurement{V1}}, ::Type{Dual{T, V2, N}}) where {V1<:AbstractFloat, T, V2, N}
    Vx = promote_rule(Measurement{V1},V2)
    return Dual{T , Vx, N}
end

function overload_ambiguous_binary(M,f)
    Mf = :($M.$f)
    return quote
        @inline function $Mf(x::Dual{Tx}, y::Measurement) where {Tx}
            ∂y = Dual{Tx}(y)
            $Mf(x,∂y)
        end
        
        @inline function $Mf(x::Measurement,y::Dual{Ty}) where {Ty}
            ∂x = Dual{Ty}(x)
            $Mf(∂x,y)
        end
    end
end

#use DiffRules.jl rules

for (M, f, arity) in DiffRules.diffrules(filter_modules = nothing)
    if (M, f) in ((:Base, :^), (:NaNMath, :pow))
        continue  # Skip methods which we define elsewhere.
    elseif !(isdefined(@__MODULE__, M) && isdefined(getfield(@__MODULE__, M), f))
        continue  # Skip rules for methods not defined in the current scope
    end
    if arity == 2
        eval(overload_ambiguous_binary(M,f))
    else
        # error("ForwardDiff currently only knows how to autogenerate Dual definitions for unary and binary functions.")
        # However, the presence of N-ary rules need not cause any problems here, they can simply be ignored.
    end
end
