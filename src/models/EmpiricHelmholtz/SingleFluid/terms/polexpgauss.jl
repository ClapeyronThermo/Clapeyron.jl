struct PolExpGaussTerm  <: MultiParameterTerm
    iterators::MVector{3,UnitRange{Int}}
    n::Vector{Float64}
    t::Vector{Float64}
    d::Vector{Float64}
    l::Vector{Float64}
    g::Vector{Float64}
    eta::Vector{Float64}
    beta::Vector{Float64}
    gamma::Vector{Float64}
    epsilon::Vector{Float64}
    function PolExpGaussTerm(n,t,d,l = Float64[],g = ones(length(l)),
        eta = Float64[],beta = Float64[],gamma = Float64[], epsilon = Float64[])
        iterators = MVector((1:0,1:0,1:0))
        param = new(iterators,n,t,d,l,g,eta,beta,gamma,epsilon)
        _calc_iterators!(param)
        return param
    end
end

PolExpGaussTerm() = PolExpGaussTerm(Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[])

active_term(term::PolExpGaussTerm) = !iszero(length(term.n))

function Base.empty!(term::PolExpGaussTerm)
    Base.empty!(term.n)
    Base.empty!(term.t)
    Base.empty!(term.d)
    Base.empty!(term.l)
    Base.empty!(term.g)
    Base.empty!(term.eta)
    Base.empty!(term.beta)
    Base.empty!(term.gamma)
    Base.empty!(term.epsilon)
    _calc_iterators!(term)
    return term
end

function a_term(term::PolExpGaussTerm,δ,τ,lnδ,lnτ,_0)
    αᵣ = _0
    if active_term(term)
        n,t,d = term.n,term.t,term.d
        k_pol,k_exp,k_gauss = term.iterators

        #strategy for storing.
        #n, t, d, gauss values, always require views
        #l, b does not require views. they are used just once.

        #Polynomial terms
        n_pol = view(n,k_pol)
        t_pol = view(t,k_pol)
        d_pol = view(d,k_pol)
        αᵣ += term_ar_pol(δ,τ,lnδ,lnτ,αᵣ,n_pol,t_pol,d_pol)

        #Exponential terms.
        if length(k_exp) != 0
            l,g = term.l,term.g
            n_exp = view(n,k_exp)
            t_exp = view(t,k_exp)
            d_exp = view(d,k_exp)
            αᵣ += term_ar_exp(δ,τ,lnδ,lnτ,αᵣ,n_exp,t_exp,d_exp,l,g)
        end

        #Gaussian bell-shaped terms
        η,β,γ,ε = term.eta,term.beta,term.gamma,term.epsilon
        if length(k_gauss) != 0
            n_gauss = view(n,k_gauss)
            t_gauss = view(t,k_gauss)
            d_gauss = view(d,k_gauss)
            αᵣ += term_ar_gauss(δ,τ,lnδ,lnτ,αᵣ,n_gauss,t_gauss,d_gauss,η,β,γ,ε)
        end
    end
    return αᵣ
end