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
        #note: the eta term in the paper has the opposite sign that the one used for the parser.
        αᵣ += n[k]*exp(lnδ*d[k] + lnτ*t[k] + η[k]*Δδ*Δδ  + 1/(β[k]*Δτ*Δτ + b[k]))
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
        n[k]*δ*Ψ*Δ^b[k]
        αᵣ += n[k]*δ*Ψ*Δ^b[k]
    end
    return αᵣ
end

@inline function _f0_gpe(τ,lnτ,_0,n,t,c,d)
    αᵣ = zero(_0)
    for k in eachindex(n)
        αᵣ += n[k]*log(c[k] + d[k]*exp(t[k]*τ))
    end
    return αᵣ
end

@inline function _f0_power(τ,logτ,_0,n,t)
    αᵣ = zero(_0)
    for k in eachindex(n)
        α₀ += n[k]*exp(logτ*t[k])
    end
    return αᵣ
end

#=
function iapws95_f0(δ,τ)
    nδ1,nδ2 = (-0.14874640856724*δ, 0.31806110878444*δ)
    _τ = (τ-1.0)^2
    _δ = (δ-1.0)^2
    Θ = (1.0-τ) + 0.32*_δ^(1.6666666666666667)
    Δ = Θ^2 + 0.2*_δ^3.5
    Ψ1 = exp(- 28*_δ - 700*_τ)
    Ψ2 = exp(- 32*_δ - 800*_τ)
    Δb1 = Δ^0.85
    Δb2 = Δ^0.95
    @show Δb1,Ψ1
    @show Δb2,Ψ2
    res = nδ1*Δb1*Ψ1 + nδ2*Δb2*Ψ2
    @show res
    return res
end=#