function lb_volume(model::SAFTModel, z = SA[1.0])
    m = model.params.segment.values
    σ = model.params.sigma.values
    m_idx = linearidx(m)
    σ_idx = linearidx(σ)
    val = zero(Base.promote_eltype(m,σ,z))
    for i in 1:length(model)
        mi = m[m_idx[i]]
        σi = σ[σ_idx[i]]
        val += z[i]*mi*σi*σi*σi
    end
    return π/6*N_A*val
end

"""
    ck_diameter(model, T, z, k1 = 0.12, k2 = 3.0)

Chen and Kregleswski efective diameter.
```
dᵢ = σᵢ*(1 - k1*exp(-k2ᵢ/ T))
```
"""
function ck_diameter(model, T, z, k1 = 0.12, k2 = 3.0)
    ϵ = model.params.epsilon.values
    σ = model.params.sigma.values
    di = zeros(Base.promote_eltype(model,T,z,k1,k2),length(model))
    for i in eachindex(di)
        di[i] = σ[i,i]*(1 - k1*exp(-k2*ϵ[i,i]/ T))
    end
    return di
end

function ck_diameter(model, T, z::SingleComp,k1 = 0.12, k2 = 3.0)
    ϵ = only(model.params.epsilon.values)
    σ = only(model.params.sigma.values)
    return SA[σ*(1 - k1*exp(-k2*ϵ/T))]
end

function ζ0123(model, V, T, z, _d=@f(d),m = model.params.segment.values)
    #N_A*π/6/V * sum(z[i]*m[i]*@f(d,i)^n for i ∈ @comps)
    _0 = zero(V+T+first(z)+one(eltype(model)))
    d_idx = linearidx(_d)
    m_idx = linearidx(m)
    ζ0,ζ1,ζ2,ζ3 = _0,_0,_0,_0
    @inbounds for i ∈ 1:length(z)
        di =_d[d_idx[i]]
        xS = z[i]*m[m_idx[i]]
        ζ0 += xS
        ζ1 += xS*di
        ζ2 += xS*di*di
        ζ3 += xS*di*di*di
    end
    c = π/6*N_A/V
    ζ0,ζ1,ζ2,ζ3 = c*ζ0,c*ζ1,c*ζ2,c*ζ3
    return ζ0,ζ1,ζ2,ζ3
end

function ζ(model, V, T, z, n, _d = @f(d),m = model.params.segment.values)
    #N_A*π/6/V * sum(z[i]*m[i]*@f(d,i)^n for i ∈ @comps)
    _0 = zero(V+T+first(z)+one(eltype(model)))
    d_idx = linearidx(_d)
    m_idx = linearidx(m)
    ζn = _0
    @inbounds for i ∈ 1:length(z)
        di =_d[d_idx[i]]
        xS = z[i]*m[m_idx[i]]
        ζn += xS*di^n
    end
    c = π/6*N_A/V
    ζn = c*ζn
    return ζn
end

g_hs_ij(d, ζ2, ζ3, i::Integer, j::Integer) = g_hs_ij(d[i], d[j], ζ2, ζ3)

function g_hs_ij(di, dj, ζ2, ζ3)
    dij = di*dj/(di+dj)
    ζ3inv = 1/(1-ζ3)
    dζ = dij*ζ3inv
    return ζ3inv + dζ*3ζ2*ζ3inv + 2*(dζ*dζ*ζ2*ζ2)*ζ3inv
end

#=
when you evaluate an EoS at zero volume, the hard sphere term diverges.
teqp implements an alternate expression to evaluate the boublik-mansoori-carnahan-starling
hard-sphere term such as the derivatives are correct

=#
function bmcs_hs_zero_v(model,V,T,z,_d = @f(d),m = model.params.segment.values)
    _0 = zero(V+T+first(z)+one(eltype(model)))
    d_idx = linearidx(_d)
    m_idx = linearidx(m)
    ζ0V,ζ1V,ζ2V,ζ3V = _0,_0,_0,_0
    @inbounds for i ∈ 1:length(z)
        di =_d[d_idx[i]]
        xS = z[i]*m[m_idx[i]]
        ζ0V += xS
        ζ1V += xS*di
        ζ2V += xS*di*di
        ζ3V += xS*di*di*di
    end
    c = π/6*N_A
    D0,D1,D2,D3 = c*ζ0V,c*ζ1V,c*ζ2V,c*ζ3V
    ρ = 1/V
    ζ2,ζ3 = ρ*ζ2V,ρ*ζ3V
    Δζ3 = 1.0 - ζ3
    logΔζ3 = log(Δζ3)
    return  3.0*D1/D0*ζ2/Δζ3
            + D2*D2*ζ2/(D3*D0*Δζ3*Δζ3)
            - logΔζ3
            + (D2*D2*D2)/(D3*D3*D0)*logΔζ3
end

"""
    packing_fraction(model, V, T, z)
    packing_fraction(model,data)

Calculates the packing fraction, defined as:
```
π/6*N_A/v * ∑xᵢmᵢdᵢ^3
```
"""
function packing_fraction(model, V, T, z)
    return ζ(model,V,T,z,3)
end

function packing_fraction(model, V, T, z, _d, m)
    return ζ(model,V,T,z,3,_d,m)
end

#fast getter in case you already calculated the packing fraction.
#overload in the following way:
# packing_fraction(model::MyModel,data::Tuple)
# packing_fraction(model,data) = nothing

function x0_crit_pure(model::SAFTModel)
    lb_v = lb_volume(model)
    (2.0, log10(lb_v/0.3))
end

function saft_lorentz_berthelot(params)
    k = get(params,"k",nothing)
    l = get(params,"l",nothing)
    sigma,epsilon = params["sigma"],params["epsilon"]
    params["sigma"] = sigma_LorentzBerthelot(sigma, l)
    params["epsilon"] = epsilon_LorentzBerthelot(epsilon, k)
    return params
end

function T_scale(model::SAFTModel,z)
    ϵ = model.params.epsilon.values
    return prod(ϵ[i,i]^z[i] for i in 1:length(z))^(1/sum(z))
end

function T_scales(model::SAFTModel)
    ϵ =diagvalues(model.params.epsilon)
end

function p_scale(model::SAFTModel,z)
    ϵ = model.params.epsilon.values
    σ = model.params.sigma.values
    V = zero(Base.promote_eltype(ϵ,σ,z))
    T = zero(Base.promote_eltype(ϵ,σ,z))
    for i in 1:length(z)
        zi = z[i]
        V += zi*N_A*σ[i,i]^3
        T += zi*ϵ[i,i]
    end
    return Rgas(model)*T/V
end

function antoine_coef(model::SAFTModel)
    m = model.params.segment.values[1]
    A = 2.3461144513376593+0.27679968565666935*m
    B = exp(1.7330494260220226 + 0.6185684341246401*log(m))
    C = 0.018524160155803788 - 0.19222021003570597*log(m)
    return A,B,C
end

#recombine! utilities
function recombine_saft!(model::SAFTModel,k = nothing,l = nothing)
    sigma = model.params.sigma
    epsilon = model.params.epsilon
    sigma = sigma_LorentzBerthelot!(sigma,l)
    epsilon = epsilon_LorentzBerthelot!(epsilon,k)
    recombine_assoc!(model)
    return model
end

