struct TillnerRothFriendDeparture <: MultiFluidDepartureModel
    components::Vector{String}
    nh3_idx::SpecialComp
end

function multiparameter_a_res(model::MultiFluid,V,T,z,departure::TillnerRothFriendDeparture,δ,τ,∑z = sum(z))
    nh3 = departure.nh3_idx
    lnδ = log(δ)
    lnτ = log(τ)
    aᵣ = multiparameter_a_res0(model,V,T,z,δ,τ,lnδ,lnτ,∑z)
    _0 = zero(aᵣ)
    isone(length(z)) && return aᵣ
    ∑z⁻¹ = 1/∑z
    #a = (-1.855822E-02,5.258010E-02,3.552874E-10,5.451379E-06,-5.998546E-13,-3.687808E-06,0.2586192,-1.368072E-08,1.226146E-02,-7.181443E-02,9.970849E-02,1.0584086E-03,-0.1963687,-0.7777897)
    #t = (1.5,0.5,6.5,1.75,15.0,6.0,-1.0,4.0,3.5,0.0,-1.0,8.0,7.5,4.0)
    #d = (4,5,15,12,12,15,4,15,4,5,6,10,6,2)
    #e = (0,1,1,1,1,2,1,1,1,1,2,2,2,2)

    a0,t0,d0 = (-0.01855822,),(1.5,),(4,)

    Δa0 = term_ar_pol(δ,τ,lnδ,lnτ,_0,a0,t0,d0)

    a1 = (0.0525801, 3.552874e-10, 5.451379e-6, -5.998546e-13, -3.687808e-6)
    t1 = (0.5, 6.5, 1.75, 15.0, 6.0)
    d1 = (5, 15, 12, 12, 15)
    e1 = (1, 1, 1, 1, 2)
    g1 = (1, 1, 1, 1, 1)
    Δa1 = term_ar_exp(δ,τ,lnδ,lnτ,_0,a1,t1,d1,e1,g1)

    a2 = (0.2586192, -1.368072e-8, 0.01226146, -0.07181443, 0.09970849, 0.0010584086, -0.1963687)
    t2 = (-1.0, 4.0, 3.5, 0.0, -1.0, 8.0, 7.5)
    d2 = (4, 15, 4, 5, 6, 10, 6)
    e2 = (1, 1, 1, 1, 2, 2, 2)
    g2 = (1, 1, 1, 1, 1, 1, 1)
    Δa2 = term_ar_exp(δ,τ,lnδ,lnτ,_0,a2,t2,d2,e2,g2)

    a3,t3,d3,e3 = (-0.7777897, 4.0, 2, 2)
    Δa3 = term_ar_exp(δ,τ,lnδ,lnτ,_0,a3,t3,d3,e3,1)

    γ = 0.5248379
    if departure.nh3_idx[] == 1
        x,xw = z[1]*∑z⁻¹,z[2]*∑z⁻¹
    else
        x,xw = z[2]*∑z⁻¹,z[1]*∑z⁻¹
    end
    Δa = (1-xw)*(1-x^γ)*(Δa0 + Δa1 + x*Δa2 + x*x*Δa3)
    return aᵣ + Δa
end

export TillnerRothFriendDeparture