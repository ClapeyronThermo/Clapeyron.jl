@inline function term_ar_pol(δ,τ,lnδ,lnτ,_0,n,t,d)
    αᵣ = zero(_0)
    for k in eachindex(n)
        αᵣ += n[k]*exp(lnδ*d[k] + lnτ*t[k])
    end
    return αᵣ
end

function term_ar_exp(δ,τ,lnδ,lnτ,_0,n,t,d,l,g)
    αᵣ = zero(_0) 
    for k in eachindex(n)
        dpart = lnδ*d[k] - g[k]*δ^l[k]
        αᵣ += n[k]*exp(dpart + lnτ*t[k])
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

#we use the specific GERG term instead of transforming into power terms because
#`logcosh` and `logabssinh` are one eval,vs 4 evals resulting in the transformation.
@inline function term_a0_gerg2008(τ,logτ,_0,n,v)
    α₀ = zero(_0)
    @inbounds for i in 2:2:length(n)
        nᵢ, ϑᵢ = n[i], v[i] #even terms
        nⱼ, ϑⱼ = n[i-1], v[i-1] #odd terms
        iszero(nᵢ) || (α₀ -= nᵢ*EoSFunctions.logcosh(ϑᵢ*τ))
        iszero(nⱼ) || (α₀ += nⱼ*EoSFunctions.logabssinh(ϑⱼ*τ))
    end
    return α₀
end

function _Cp0_constant_parse(c,Tc,T0)
    #c - cT0/Tc*τ + c*(log(τ/τ0))
    #c - c*τ0inv*τ + c*log(τ*τ0inv)
    #c - c*τ0inv*τ + c*log(τ) + c*(τ0inv)
    τ0inv = T0/Tc
    a1 = c*(1 + log(τ0inv))
    a2 = -c*τ0inv
    c0 = c
    return a1,a2,c0
end

function _Cpi_power_parse(c,t,Tc,T0)
    a1 = c*(T0^t)/t
    a2 = -c * (T0^(t+1)) / (Tc * (t + 1))
    ni = -c * (Tc^t) / (t * (t + 1))
    ti = -t
    return a1,a2,ni,ti
end