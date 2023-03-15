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
        αᵣ += n[k]*exp(lnδ*d[k] + lnτ*t[k] - η[k]*(δ-ε[k])^2  - β[k]*(τ-γ[k])^2)
    end
    return αᵣ
end

@inline function _fr1_gerg2008(δ,τ,lnδ,lnτ,_0,n,t,d,η,β,γ,ε)
    αᵣ = zero(_0)
    for k in eachindex(n)
        αᵣ += n[k]*exp(lnδ*d[k] + lnτ*t[k] - η[k]*(δ-ε[k])^2  - β[k]*(δ-γ[k]))
    end
    return αᵣ
end

@inline function _fr1_gao(δ,τ,lnδ,lnτ,_0,n,t,d,η,β,γ,ε,b)
    αᵣ = zero(_0)
    for k in eachindex(n)
        αᵣ += n[k]*exp(lnδ*d[k] + lnτ*t[k] - η[k]*(δ-ε[k])^2  + 1/(β[k]*(τ-γ[k])^2 + b[k]))
    end
    return αᵣ
end