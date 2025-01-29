using Clapeyron: viewlast, viewn, viewfirst
function positive_linesearch(v, δ, α0 = 1.0 ; τ = 1.0, decay = 0.5, tol = 1e-10, s = 1.0)
    done = false
    l = length(v)
    α = α0*oneunit(Base.promote_eltype(v, δ, τ, s))

    iszero(l) && return α
        for i in 1:length(v)
            vi,δi = v[i],δ[i]
            if vi + α * s * δi < (1 - τ) * vi
                α = min(α*decay,(-τ*vi)/(s*δi))
            end
        end
    return α*decay
end

function backtracking_linesearch!(Θ,F,X,d,Θ0,Xnew,α = 1.0;tol = 1e-10, decay = 0.5, ignore = nothing)
    done = false
    Θx = Θ0
    while !done
        if ignore === nothing
            Xnew .= X .+ α .* d
        else
            for i in 1:length(X)
                if !ignore[i]
                    Xnew[i] = X[i] + α * d[i]
                end
            end
        end
        Θx = Θ(F,Xnew)
        if Θx <= Θ0
            done = true
        else
            if α < tol || !isfinite(Θx)
                return zero(α)/zero(α),Θx
            end
            α *= decay
        end
    end
    return α,Θx
end

function remove_slacks!(F,J,slacks::AbstractVector{Bool})
    for i in 1:length(slacks)
        if slacks[i]
            remove_slacks!(F,J,i)
        end
    end
end

function remove_slacks!(F,J,i::Int)
    if F !== nothing
        F[i] = 0
    end

    if J !== nothing
        J1 = @view(J[i,:])
        J2 = @view(J[:,i])
        J1 .= 0
        J2 .= 0
        J[i,i] = 1
    end
end
#=
"""
    npipm(Λ!,G!,H!,x0,l0,g0,h0))

Non-Parametric Interior point method, used to solve systems of the form:
```
l = Λ(x) = 0
g = G!(x) = 0
h = H!(x) = 0
g .* h = 0
g >= 0
h >= 0
```
"""
function npipm(Λ!::L,G!::G,H!::H,x,m::Int,positives = 1:0;rtol = 1e-12,atol = 1e-10,max_iters = 50) where {L,G,H}
    l = length(x)
    lx = l + 2*m + 1
    𝕏 = similar(x,lx)
    𝕏[1:length(x)] .= x
    ∇𝔽 = similar(x,lx,lx)
    𝔽 = similar(𝕏)
    𝕕 = similar(𝕏)
    𝕏old = similar(𝕏)
    𝕏old .= 0
    vwv = viewlast(𝔽,3*m+1)
    x_vwv = viewlast(𝕏,3*m+1)
    V0 = viewn(vwv,m,1)
    W0 = viewn(vwv,m,2)

    G!(V0,x)
    H!(W0,x)
    @show V0
    XV0 = viewn(x_vwv,m,1)
    XW0 = viewn(x_vwv,m,2)

    XW0 .= W0 .+ 1
    XV0 .= V0 .+ 1
    @show XW0,XV0
    η = 0.5
    u = 1.0
    v0 = npipm_v0(XV0,XW0,u,m,η)
    @show v0
    𝕏[end] = v0
    𝔽!(F,X) = npipm_neq(Λ!,G!,H!,X,F,m,u,η)
    Θ(f) = 0.5*dot(f,f)
    Θ(_f,z) = Θ(𝔽!(_f,z))
    Θx = Θ(𝔽,𝕏)
    @show Θx
    𝔽norm = sqrt(2*Θx)
    𝔽norm_old = Inf*𝔽norm
    config = ForwardDiff.JacobianConfig(𝔽!,𝔽,𝕏)
    piv = zeros(Int,lx)
    converged = false
    nan_converged = false
    for i in 1:50
        converged && break
        nan_converged && break
        ForwardDiff.jacobian!(∇𝔽,𝔽!,𝔽,𝕏,config)
        display(∇𝔽)
        for i in 1:length(𝔽)
            if abs(𝔽[i]) < eps()
                𝔽[i] = 0.0
            end
        end
        #display(∇𝔽)
        #lu = Solvers.unsafe_LU!(∇𝔽,piv)
        𝕕 .= -𝔽
        #ldiv!(lu,𝕕) #s .= J\F
        𝕕 .= sparse(∇𝔽)\𝔽

        𝕏old .= 𝕏
        @show 𝕏old
        @show 𝕕
        #backtrack positive variables in the original X
        #α1 = positive_linesearch(@view(𝕏[positives]),@view(𝕕[positives]))
        α1 = 1.0
        #backtrack V,W,v to positive
        α2 = positive_linesearch(viewlast(𝕏,2*m+1),viewlast(𝕕,2*m+1),α1)

        #backtrack acording to 0.5*||F||^2
        α3,Θx = backtracking_linesearch!(Θ,𝔽,𝕏old,𝕕,Θx,𝕏,α2)
        @show α1,α2,α3 
        𝔽norm_old = 𝔽norm
        𝔽norm = sqrt(2*Θx)
        #there could be the case that some fugacity coefficients cause problems
        #the algorithm converges, but the norm does not.
        
        #convergence criteria
        abs(𝔽norm_old-𝔽norm) < 1e-8         && (converged = true)
        norm(𝔽,Inf) < rtol                  && (converged = true)
        Solvers.dnorm(𝕏,𝕏old,Inf) < atol    && (converged = true)

        nan_converged = !all(isfinite,𝕏)
    end
    x .= @view 𝕏[1:l]
    return x
end

function npipm_neq(Λ!::L,G!::G,H!::H,𝕏,𝔽,m,u,η) where {L,G,H}
    𝔽 .= 0
    X =  @view 𝕏[1:end - 2m - 1]
    F = @view 𝔽[1:end - 2m - 1]
    Λ!(F,X)
    vwv,F_vwv = viewlast(𝕏,3*m+1),viewlast(𝔽,3*m+1)
    GX,HX = viewn(F_vwv,m,1),viewn(F_vwv,m,2)
    V,W = viewn(vwv,m,1),viewn(vwv,m,2)
    Z = viewn(F_vwv,m,3)
    v = last(𝕏)
    G!(GX,X)
    H!(HX,X)
    GX .-= V #G(X) - V
    HX .-= W #H(X) - W
    Z .= V .* W .- v #  V ⊙ W .- v
    f = npipm_fv(V,W,u,m,η,v)
    𝔽[end] = f
    return 𝔽
end

function norm_minus2(W)
    w = zero(eltype(W))
    _0 = zero(eltype(W))
    @inbounds for i in eachindex(W)
        w += abs2(min(W[i],_0))
    end
    return w
end

function npipm_v0(V,W,u,m,η)
    V2 = norm_minus2(V)
    W2 = norm_minus2(W)
    VW = dot(V,W)
    return VW/m
    VW⁺ = max(VW,zero(VW))
    VW2 = VW⁺*VW⁺
    k = 0.5*(V2 + W2) + 0.5*(u/(m*m))*VW2
    @show k
    v = 0.5*(sqrt(η*η - 4*k) - η)
end

function npipm_fv(V,W,u,m,η,v)
    V2 = norm_minus2(V)
    W2 = norm_minus2(W)
    VW = dot(V,W)
    VW⁺ = max(VW,zero(VW))
    VW2 = VW⁺*VW⁺
    f = 0.5*(V2 + W2) + 0.5*(u/(m*m))*VW2 + η*v + v*v
end
=#