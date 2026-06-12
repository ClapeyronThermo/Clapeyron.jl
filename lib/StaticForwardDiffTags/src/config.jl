function ForwardDiff.DerivativeConfig(f::F,
                          y::AbstractArray{Y},
                          x::X,
                          tag::T = maketag(f, X)) where {F<:SDiffFunction,X<:Real,Y<:Real,T}
    return ForwardDiff.DerivativeConfig(f,y,x,tag)
end

function ForwardDiff.GradientConfig(f::F,
                        x::AbstractArray{V},
                        chunk::Chunk{N} = Chunk(x),
                        tag::T = maketag(f, V)) where {F<:SDiffFunction,V,N,T}
    return ForwardDiff.GradientConfig(f,x,chunk,tag)
end

function JacobianConfig(f::F,
                        x::AbstractArray{V},
                        chunk::Chunk{N} = Chunk(x),
                        tag::T = maketag(f, V)) where {F<:SDiffFunction,V,N,T}
    return ForwardDiff.JacobianConfig(f,x,chunk,tag)
end

function ForwardDiff.JacobianConfig(f::F,
                        y::AbstractArray{Y},
                        x::AbstractArray{X},
                        chunk::Chunk{N} = Chunk(x),
                        tag::T = maketag(f, X)) where {F<:SDiffFunction,Y,X,N,T}
    return ForwardDiff.JacobianConfig(f,y,x,chunk,tag)
end

function ForwardDiff.HessianConfig(f::F,
                       x::AbstractArray{V},
                       chunk::Chunk = Chunk(x),
                       tag = maketag(f, V)) where {F<:SDiffFunction,V}
    return ForwardDiff.HessianConfig(f,x,chunk,tag)
end

function ForwardDiff.HessianConfig(f::F,
                       result::DiffResults.DiffResult,
                       x::AbstractArray{V},
                       chunk::Chunk = Chunk(x),
                       tag = Tag(f, V)) where {F<:SDiffFunction,V}
    ForwardDiff.HessianConfig(f,result,x,chunk,tag)
end

