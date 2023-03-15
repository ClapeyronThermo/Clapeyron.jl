#Always try to minimize ^ operations, those are really really costly.
@inline function _fr1_pol(δ,τ,lnδ,lnτ,_0,n,t,d)
    αᵣ = zero(_0)
    for k in eachindex(n)
        αᵣ += n[k]*exp(lnδ*d[k] + lnτ*t[k])
    end
    return αᵣ
end

@inline function _fr1_exp(δ,τ,lnδ,lnτ,_0,n,t,d,l)
    αᵣ = zero(_0)
    for k in eachindex(n)
        αᵣ += n[k]*exp(lnδ*d[k] + lnτ*t[k] - δ^l[k])
    end
    return αᵣ
end

@inline function _fr1_gauss(δ,τ,lnδ,lnτ,_0,n,t,d,η,β,γ,ε)
    αᵣ = zero(_0)
    for k in eachindex(n)
        Δδ = δ-ε[k]
        Δτ = τ-γ[k]
        αᵣ += n[k]*exp(lnδ*d[k] + lnτ*t[k] - η[k]*Δδ*Δδ  - β[k]*Δτ*Δτ)
    end
    return αᵣ
end

@inline function _fr1_gerg2008(δ,τ,lnδ,lnτ,_0,n,t,d,η,β,γ,ε)
    αᵣ = zero(_0)
    for k in eachindex(n)
        Δδ = δ-ε[k]
        αᵣ += n[k]*exp(lnδ*d[k] + lnτ*t[k] - η[k]*Δδ*Δδ  - β[k]*Δδ)
    end
    return αᵣ
end

@inline function _fr1_gao(δ,τ,lnδ,lnτ,_0,n,t,d,η,β,γ,ε,b)
    αᵣ = zero(_0)
    for k in eachindex(n)
        Δδ = δ-ε[k]
        Δτ = τ-γ[k]
        αᵣ += n[k]*exp(lnδ*d[k] + lnτ*t[k] - η[k]*Δδ*Δδ  + 1/(β[k]*Δτ*Δτ + b[k]))
    end
    return αᵣ
end

@inline function _fr1_na(δ,τ,lnδ,lnτ,_0,A,B,C,D,a,b,β,n)
    αᵣ = zero(_0)
    Δδ = δ-1
    Δτ = τ-1
    logΔδ2 = log(Δδ*Δδ)
    for k in eachindex(n)
        Ψ = exp(-C[k]*Δδ*Δδ - D[k]*Δτ*Δτ)
        Θ = -Δτ + A[k]*exp(logΔδ2/(2*β[k]))
        Δ = Θ*Θ + B[k]*exp(logΔδ2*a[k])
        αᵣ += n[k]*δ*Ψ*Δ^b[i]
    end
    return αᵣ
end

@inline function _f0_gpe(τ,_0,n,t,c,d)
    αᵣ = zero(_0)
    for k in eachindex(n)
        αᵣ += n[k]*exp(lnδ*d[k] + lnτ*t[k])
    end
    return αᵣ
end

@inline function _f0_power(τ,_0,n,t)
    αᵣ = zero(_0)
    logτ = log(τ)
    for i in eachindex(n)
        α₀ += n[i]*exp(logτ*t[i])
    end
    return αᵣ
end