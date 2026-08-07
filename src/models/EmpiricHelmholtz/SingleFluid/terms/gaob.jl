struct GaoBTerm <: MultiParameterTerm
    n::Vector{Float64}
    t::Vector{Float64}
    d::Vector{Float64}
    eta::Vector{Float64}
    beta::Vector{Float64}
    gamma::Vector{Float64}
    epsilon::Vector{Float64}
    b::Vector{Float64}
    function GaoBTerm(n,t,d,eta,beta,gamma,epsilon,b)
        @assert length(eta) == length(beta) == length(gamma) == length(epsilon) == length(b)
        @assert length(eta) == length(n) == length(t) == length(d)
        return new(n,t,d,eta,beta,gamma,epsilon,b)
    end
end

GaoBTerm() = GaoBTerm(Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[])

active_term(term::GaoBTerm) = !iszero(length(term.n))

function Base.empty!(term::GaoBTerm)
    Base.empty!(term.n)
    Base.empty!(term.t)
    Base.empty!(term.d)
    Base.empty!(term.eta)
    Base.empty!(term.beta)
    Base.empty!(term.gamma)
    Base.empty!(term.epsilon)
    Base.empty!(term.b)
end

function a_term(term::GaoBTerm,δ,τ,lnδ,lnτ,_0)
    if active_term(term)
        n = term.n
        t = term.t
        d = term.d
        η = term.eta
        β = term.beta
        γ = term.gamma
        ε = term.epsilon
        b = term.b
        αᵣ = term_ar_gaob(δ,τ,lnδ,lnτ,_0,n,t,d,η,β,γ,ε,b)
    else
        αᵣ = _0
    end
    return αᵣ
end