@inline function term_ar_pol(δ,τ,lnδ,lnτ,_0,n,t,d)
    αᵣ = zero(_0)
    for k in eachindex(n)
        αᵣ += n[k]*exp(lnδ*d[k] + lnτ*t[k])
    end
    return αᵣ
end

@inline function term_ar_gauss(δ,τ,lnδ,lnτ,_0,n,t,d,η,β,γ,ε)
    αᵣ = zero(_0)
    for k in eachindex(n)
        Δδ = δ-ε[k]
        Δτ = τ-γ[k]
        αᵣ += n[k]*exp(lnδ*d[k] + lnτ*t[k] - η[k]*Δδ*Δδ  - β[k]*Δτ*Δτ)
    end
    return αᵣ
end
#TODO: transform to gauss form?
#(described in EOS-LNG paper)
@inline function term_ar_gerg2008(δ,τ,lnδ,lnτ,_0,n,t,d,η,β,γ,ε)
    αᵣ = zero(_0)
    for k in eachindex(n)
        Δδ = δ-ε[k]
        αᵣ += n[k]*exp(lnδ*d[k] + lnτ*t[k] - η[k]*Δδ*Δδ  - β[k]*(δ-γ[k]))
    end
    return αᵣ
end

@inline function term_ar_gaob(δ,τ,lnδ,lnτ,_0,n,t,d,η,β,γ,ε,b)
    αᵣ = zero(lnδ+lnτ)
    for k in eachindex(n)
        Δδ = δ-ε[k]
        Δτ = τ-γ[k]
        #note: the eta term in the paper has the opposite sign that the one used for the parser.
        αᵣ += n[k]*exp(lnδ*d[k] + lnτ*t[k] + η[k]*Δδ*Δδ  + 1/(β[k]*Δτ*Δτ + b[k]))
    end
    return αᵣ
end

@inline function term_ar_na(δ,τ,lnδ,lnτ,_0,A,B,C,D,a,b,β,n)
    αᵣ = zero(_0)
    Δδ = δ-1
    Δτ = τ-1
    logΔδ2 = log(Δδ*Δδ)
    for k in eachindex(n)
        Ψ = exp(-C[k]*Δδ*Δδ - D[k]*Δτ*Δτ)
        Θ = -Δτ + A[k]*exp(logΔδ2/(2*β[k]))
        Δ = Θ*Θ + B[k]*exp(logΔδ2*a[k])
        n[k]*δ*Ψ*Δ^b[k]
        αᵣ += n[k]*δ*Ψ*Δ^b[k]
    end
    return αᵣ
end

function term_ar_assoc2b(δ,τ,lnδ,lnτ,_0,ε,κ,a,m,v̄ₙ)
    η = v̄ₙ*δ
    g = 0.5*(2 - η)/(1 - η)^3
    Δ = g*(exp(ε*τ) - 1)*κ
    X = 2 / (sqrt(1 + 4 * Δ * δ) + 1)
    return m * a * ((log(X) - X / 2.0 + 0.5))
end

function term_ar_exp(δ,τ,lnδ,lnτ,_0,n,t,d,l,g)
    αᵣ = zero(_0)
    for k in eachindex(n)
        αᵣ += n[k]*exp(lnδ*d[k] + lnτ*t[k] - g[k]*δ^l[k])
    end
    return αᵣ
end

#ideal terms

@inline function term_a0_gpe(τ,lnτ,_0,n,t,c,d)
    α₀ = zero(_0)
    for k in eachindex(n) 
        #αᵣ += n[k]*log(muladd(d[k],exp(-t[k]*τ),c[k]))
        α₀ += n[k]*log(d[k]*exp(t[k]*τ) + c[k])
    end
    return α₀
end

@inline function term_a0_power(τ,logτ,_0,n,t)
    α₀ = zero(_0)
    for k in eachindex(n)
        α₀ += n[k]*exp(logτ*t[k])
    end
    return α₀
end

@inline function term_a0_gerg2008(τ,logτ,_0,n,v)
    α₀ = zero(_0)
    n₁,n₂,n₃,n₄ = n[1],n[2],n[3],n[4]
    ϑ₁,ϑ₂,ϑ₃,ϑ₄ = v[1],v[2],v[3],v[4]
    iszero(n₁) || (α₀ += n₁*LogExpFunctions.logabssinh(ϑ₁*τ))
    iszero(n₂) || (α₀ -= n₂*LogExpFunctions.logcosh(ϑ₂*τ))
    iszero(n₃) || (α₀ += n₃*LogExpFunctions.logabssinh(ϑ₃*τ))
    iszero(n₄) || (α₀ -= n₄*LogExpFunctions.logcosh(ϑ₄*τ))
    return α₀
end
