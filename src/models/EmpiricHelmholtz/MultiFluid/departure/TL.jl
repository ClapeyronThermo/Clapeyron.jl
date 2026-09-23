struct TillnerRothFriendDeparture <: MultiFluidDepartureModel end
TillnerRothFriendDeparture(::Any, ::Any) = TillnerRothFriendDeparture()

function multiparameter_a_res(model::MultiFluid,V,T,z,departure::TillnerRothFriendDeparture,δ,τ,∑z = sum(z))
    lnδ = log(δ)
    lnτ = log(τ)
    aᵣ = multiparameter_a_res0(model,V,T,z,δ,τ,lnδ,lnτ,∑z)
    _0 = zero(aᵣ)
    isone(length(z)) && return aᵣ
    _, ia = _water_ammonia_indices(model)
    xa = z[ia]/∑z
    Δa = tlf_departure(model,τ, δ, xa, lnδ, lnτ)
    return aᵣ + Δa
end

function tlf_departure(model::MultiFluid, τ, δ, xa, lnδ,lnτ)
    _0 = zero(Base.promote_eltype(lnδ,lnτ,xa))
    iszero(xa) && return _0
    a = ( 0.0,        -1.855822e-02,  5.258010e-02,  3.552874e-10,  5.451379e-06,
        -5.998546e-13, -3.687808e-06,  0.2586192,    -1.368072e-08,  1.226146e-02,
        -7.181443e-02,  9.970849e-02,  1.0584086e-03, -0.1963687,    -0.7777897 )
    t = ( 0.0, 1.5, 0.5, 6.5, 1.75, 15.0, 6.0, -1.0, 4.0, 3.5, 0.0, -1.0, 8.0, 7.5, 4.0 )
    d = ( 0.0, 4.0, 5.0, 15.0, 12.0, 12.0, 15.0, 4.0, 15.0, 4.0, 5.0, 6.0, 10.0, 6.0, 2.0 )
    e = ( 0.0, 0.0, 1.0,  1.0,  1.0,  1.0,  2.0, 1.0, 1.0,  1.0, 1.0, 2.0,  2.0, 2.0, 2.0 )
    γ = 0.5248379
    S = _0
    _1 = one(_0)
    @inbounds for n in 1:15
        an = a[n]
        iszero(an) && continue                     # a[1] == 0, skip
        # τ^t_n
        term = iszero(t[n]) ? _1   : exp(t[n] * lnτ)
        # δ^d_n
        term *= iszero(d[n]) ? _1  : exp(d[n] * lnδ)
        # exp(-δ^e_n); e[n] == 0 → exp(-1) is a compile-time constant
        term *= iszero(e[n]) ? exp(-_1) : exp(-exp(e[n] * lnδ))
        # composition weight for n = 7..14
        if 7 <= n <= 13
            term *= xa
        elseif n == 14
            term *= xa * xa
        end
        S += an * term
    end
    return xa * (1.0 - exp(γ * log(xa))) * S
end

export TillnerRothFriendDeparture