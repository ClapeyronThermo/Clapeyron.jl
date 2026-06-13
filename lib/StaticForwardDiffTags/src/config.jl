function ForwardDiff.DerivativeConfig(f::F,
                          y::AbstractArray{Y},
                          x::X,
                          tag::T = maketag(f, X)) where {F<:WithContext,X<:Real,Y<:Real,T}
    duals = similar(y, Dual{T,Y,1})
    return ForwardDiff.DerivativeConfig{T,typeof(duals)}(duals)
end

function ForwardDiff.GradientConfig(f::F,
                        x::AbstractArray{V},
                        chunk::Chunk{N} = Chunk(x),
                        tag::T = maketag(f, V)) where {F<:WithContext,V,N,T}
    seeds = ForwardDiff.construct_seeds(Partials{N,V})
    duals = similar(x, Dual{T,V,N})
    return ForwardDiff.GradientConfig{T,V,N,typeof(duals)}(seeds, duals)
end

function JacobianConfig(f::F,
                        x::AbstractArray{V},
                        chunk::Chunk{N} = Chunk(x),
                        tag::T = maketag(f, V)) where {F<:WithContext,V,N,T}
    seeds = ForwardDiff.construct_seeds(Partials{N,V})
    duals = similar(x, Dual{T,V,N})
    return ForwardDiff.JacobianConfig{T,V,N,typeof(duals)}(seeds, duals)
end

function ForwardDiff.JacobianConfig(f::F,
                        y::AbstractArray{Y},
                        x::AbstractArray{X},
                        chunk::Chunk{N} = Chunk(x),
                        tag::T = maketag(f, X)) where {F<:WithContext,Y,X,N,T}

                        seeds = ForwardDiff.construct_seeds(Partials{N,X})
    yduals = similar(y, Dual{T,Y,N})
    xduals = similar(x, Dual{T,X,N})
    duals = (yduals, xduals)
    return ForwardDiff.JacobianConfig{T,X,N,typeof(duals)}(seeds, duals)
end

function ForwardDiff.HessianConfig(f::F,
                       x::AbstractArray{V},
                       chunk::Chunk = Chunk(x),
                       tag = maketag(f, V)) where {F<:WithContext,V}
    jacobian_config = ForwardDiff.JacobianConfig(f, x, chunk, tag)
    gradient_config = ForwardDiff.GradientConfig(f, jacobian_config.duals, chunk, tag)
    return ForwardDiff.HessianConfig(jacobian_config, gradient_config)
end

#TODO: use nested WithContext to define gradient function
function ForwardDiff.HessianConfig(f::F,
                       result::DiffResults.DiffResult,
                       x::AbstractArray{V},
                       chunk::Chunk = Chunk(x),
                       tag = maketag(f, V)) where {F<:WithContext,V}
    jacobian_config = ForwardDiff.JacobianConfig((f,ForwardDiff.gradient), DiffResults.gradient(result), x, chunk, tag)
    gradient_config = ForwardDiff.GradientConfig(f, jacobian_config.duals[2], chunk, tag)
    return HessianConfig(jacobian_config, gradient_config)
end

