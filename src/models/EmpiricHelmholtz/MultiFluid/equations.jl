
"""
    mixing_rule_asymetric(op, op_asym, x, p, A, A_asym)

returns an efficient implementation of:
` sum(A[i,j] * x[i] * x[j] * op(p[i],p[j]) * op_asym(x[i],x[j],A_asym[i,j])) for i = 1:n , j = 1:n)`
where `op(p[i],p[j]) == op(p[j],p[i])` , op_asym doesn't follow this symmetry.

""" 
function mixing_rule_asymetric(op, op_asym, x, p, A, A_asym)
    N = length(x)
    checkbounds(A, N, N)
    checkbounds(A_asym, N, N)
    @boundscheck checkbounds(p, N)
    @inbounds begin
        res1 = zero(eltype(x))
        for i = 1:N
            x[i] != 0 && begin
                res1 += p[i] * x[i]^2
                for j = 1:i - 1
                    res1 += 2*x[i]*x[j]*op(p[i], p[j])*A[i, j]*op_asym(x[i], x[j], A_asym[i, j])
                end
            end
        end
    end
    return res1
end

_gerg_asymetric_mix_rule(xi, xj, b) = b * (xi + xj) / (xi * b^2 + xj)

ith_index(pv::PackedVofV,i) = @inbounds begin (pv.p[i]):(pv.p[i+1]-1) end

@inline function ith_index(pol_i,exp_i,i) 
    @inbounds begin 
        kall_1::Int = pol_i[i]
        kall_end::Int = pol_i[i+1]#-1
        kexp_1::Int = exp_i[i]
        kexp_end::Int = exp_i[i+1]#-1
        lall = kall_end - kall_1
        lexp = kexp_end - kexp_1
        divider = lall - lexp
        kmid = kall_1+divider
        k1 = kall_1:(kmid - 1)
        k2 = kmid:(kall_end -1)
        kexp = kexp_1:(kexp_end-1)
        return k1,k2,kexp
    end
end

function _T_scale(model::MultiFluidModel,z=SA[1.],Σz = sum(z))
    Tc = model.properties.Tc.values
    #isone(length(z)) && return only(Tc) 
    return mixing_rule_asymetric(
        (a,b)->sqrt(a * b),
        _gerg_asymetric_mix_rule,
        z,
        Tc,
        model.ideal.gamma_T.values,
        model.ideal.beta_T.values,
    )/(Σz*Σz)
end

T_scale(model::MultiFluidModel,z=SA[1.]) = _T_scale(model,z)

p_scale(model::MultiFluidModel,z=SA[1.]) = dot(z,model.properties.pc.values)/sum(z)
 
function T_scales(model::MultiFluidModel,z=SA[1.])
    return model.properties.Tc.values
end

function _v_scale(model::MultiFluidModel,z=SA[1.],Σz = sum(z))
    vc = model.properties.vc.values
    #isone(length(z)) && return only(vc) 
    res = mixing_rule_asymetric(
        (a,b) -> ((cbrt(a) + cbrt(b))*0.5)^3,
        _gerg_asymetric_mix_rule,
        z,
        vc,
        model.ideal.gamma_v.values,
        model.ideal.beta_v.values,
    )
    return res/(Σz*Σz)
end

function lb_volume(model::MultiFluidModel,z=SA[1.])
    return dot(z,model.properties.lb_v.values)
end

function _delta(model::MultiFluidModel, rho, T, z=SA[1.],Σz = sum(z))
    vcmix = _v_scale(model,z,Σz)
    return rho * vcmix
end

function _tau(model::MultiFluidModel, rho, T, z=SA[1.],Σz = sum(z))
    Tcmix  = _T_scale(model,z,Σz)
    return Tcmix / T
end

function _f0(model::MultiFluidModel, ρ, T, z=SA[1.], Σz = sum(z))
    RR = 8.314472 / 8.314510
    #common_type = promote_type(typeof(ρ), typeof(T), eltype(x))
    _0 = zero(ρ + T + first(z))
    lnΣz = log(Σz)
    res = _0
    vc = model.properties.vc.values
    Tc = model.properties.Tc.values
    ϑ = model.ideal.theta.values
    n0 = model.ideal.n0.values
    Tinv = one(T)/T
    @inbounds for i ∈ @comps
        δi = ρ * vc[i]
        τi = Tc[i] * Tinv
        n1,n2,n3,n4,n5,n6,n7 = n0[i]
        ϑ1,ϑ2,ϑ3,ϑ4 = ϑ[i]  
        ai = n1 + n2*τi + n3*log(τi)
        iszero(n4) || (ai += n4*log(abs(sinh(ϑ1*τi))))
        iszero(n5) || (ai -= n5*log(cosh(ϑ2*τi)))
        iszero(n6) || (ai += n6*log(abs(sinh(ϑ3*τi)))) 
        iszero(n7) || (ai -= n7*log(cosh(ϑ4*τi)))
        res += z[i]*(RR*ai + log(δi) + log(z[i]) - lnΣz)
    end
    return res
end

function _fr1(model::MultiFluidModel,δ,τ,z)
    _0 = zero(promote_type(typeof(δ), typeof(τ), eltype(z)))
    res = _0
    n  = model.single.n.values
    c = model.single.c.values
    k_all = n.p
    k_exp = c.p
    nᵢ = n.v
    dᵢ = model.single.d.values.v
    tᵢ = model.single.t.values.v
    cᵢ = c.v
    lnδ = log(δ)
    lnτ = log(τ)
    @inbounds for i ∈ @comps
        ai = _0
        k1,k2,kexp = ith_index(k_all,k_exp,i)
        for k ∈ k1
            ai += nᵢ[k]*exp(lnδ*dᵢ[k] + lnτ*tᵢ[k])
            #ai += nᵢ[k]*(δ^dᵢ[k])*(τ^tᵢ[k])
        end

        for (k,k_) ∈ zip(k2,kexp)

            ai += nᵢ[k]*exp(lnδ*dᵢ[k] + lnτ*tᵢ[k]-δ^cᵢ[k_])
            #ai += nᵢ[k]*(δ^dᵢ[k])*(τ^tᵢ[k])*exp(-δ^cᵢ[k_])
        end  
        res += z[i]*ai 
    end
    return res
end

function _fr2(model::MultiFluidModel,δ,τ,z)
    _0 = zero(promote_type(typeof(δ), typeof(τ), eltype(z)))
    isone(length(z)) && return _0
    F = model.pair.Fij.values
    iszero(nnz(F)) && return _0
    res = _0
    Fij = nonzeros(F)
    n = model.pair.nij.values.storage
    η = model.pair.eta_ij.values.storage
    rows = rowvals(F)
    k_all = n.p
    k_exp = η.p
    nᵢⱼ = n.v
    tᵢⱼ = model.pair.tij.values.storage.v
    dᵢⱼ = model.pair.dij.values.storage.v
    ηᵢⱼ = η.v
    εᵢⱼ = model.pair.epsilon_ij.values.storage.v
    βᵢⱼ = model.pair.beta_ij.values.storage.v
    γᵢⱼ = model.pair.gamma_ij.values.storage.v
    lnδ = log(δ)
    lnτ = log(τ)
    @inbounds for j ∈ @comps
        for ii ∈ nzrange(F, j)
            i = rows[ii]
            Fᵢⱼ= Fij[ii]
            aij = _0
            k1,k2,kexp = ith_index(k_all,k_exp,ii)
            for k ∈ k1
                aij += nᵢⱼ[k]*exp(lnδ*dᵢⱼ[k] + lnτ*tᵢⱼ[k])
                #aij += nᵢⱼ[k]*(δ^(dᵢⱼ[k]))*(τ^(tᵢⱼ[k]))
            end
            for (k,k_) ∈ zip(k2,kexp)
                aij += nᵢⱼ[k]*(δ^(dᵢⱼ[k]))*(τ^(tᵢⱼ[k]))*
                exp(-ηᵢⱼ[k_]*(δ - εᵢⱼ[k_])^2 - βᵢⱼ[k_]*(δ -γᵢⱼ[k_]))
                #aij += nᵢⱼ[k]*exp(lnδ*dᵢⱼ[k] + lnτ*tᵢⱼ[k]
                #-ηᵢⱼ[k_]*(δ - εᵢⱼ[k_])^2 - βᵢⱼ[k_]*(δ -γᵢⱼ[k_]))
            end
           res +=z[i]*z[j]*Fᵢⱼ*aij
        end
     end
    return res
end

function a_ideal(model::MultiFluidModel, V, T, z=SA[1.0])
    Σz = sum(z)
    ρ = Σz*1.0e-3/V
    _f0(model, ρ, T, z,Σz)/Σz
end

function a_res(model::MultiFluidModel, V, T, z=SA[1.0])
    Σz = sum(z)
    invn = 1/Σz
    invn2 = invn*invn
    ρ = Σz*1.0e-3/V
    δ = _delta(model, ρ, T, z,Σz)
    τ = _tau(model, ρ, T, z,Σz)
    return (_fr1(model,δ,τ,z)*invn + _fr2(model,δ,τ,z)*invn2)
end

function eos(model::MultiFluidModel, V, T, z=SA[1.0])
    Σz = sum(z)
    invn = 1/Σz
    ρ = Σz*1.0e-3/V
    δ = _delta(model, ρ, T, z,Σz)
    τ = _tau(model, ρ, T, z,Σz)
    return  R̄*T*(_f0(model,ρ,T,z,Σz) + _fr1(model,δ,τ,z)+_fr2(model,δ,τ,z)*invn)
end

function eos_res(model::MultiFluidModel, V, T, z=SA[1.0])
    Σz = sum(z)
    invn = 1/Σz
    ρ = Σz*1.0e-3/V
    δ = _delta(model, ρ, T, z,Σz)
    τ = _tau(model, ρ, T, z,Σz)
    return  R̄*T*(_fr1(model,δ,τ,z)+_fr2(model,δ,τ,z)*invn)
end

function x0_sat_pure(model::MultiFluidModel,T)
    Ts = T_scale(model)
    vs = _v_scale(model)*0.001 #remember, vc constants in L/mol
    h = vs*5000.0
    T0 = 369.89*T/Ts
    vl = (1.0/_propaneref_rholsat(T0))*h
    vv = (1.0/_propaneref_rhovsat(T0))*h
    return (vl,vv) 
end

#Corresponding States
function x0_psat(model::MultiFluidModel,T,crit=nothing)
    Ts = T_scale(model)
    T0 = 369.89*T/Ts
    Ps = p_scale(model)
    return Ps*_propaneref_psat(T0)/4.2512e6
end

function x0_volume_liquid(model::MultiFluidModel,T,z)
    return 1.01*lb_volume(model,z)
end

function x0_volume_gas(model::MultiFluidModel,p,T,z=SA[1.])
    V = sum(z)*R̄*T/p
    return V
end

molecular_weight(model::MultiFluidModel,z=SA[1.0]) = comp_molecular_weight(mw(model),z)

mw(model::MultiFluidModel) = model.properties.Mw.values

function x0_crit_pure(model::MultiFluidModel)
    return (1.,log10(_v_scale(model)*0.001))
end

function x0_saturation_temperature(model::MultiFluidModel,p)
    p0 = p/p_scale(model)*4.2512e6
    T0 = _propaneref_tsat(p0)
    Ts = T_scale(model)
    vs = _v_scale(model)*0.001 #remember, vc constants in L/mol
    h = vs*5000.0
    vl = (1.0/_propaneref_rholsat(T0))*h
    vv = (1.0/_propaneref_rhovsat(T0))*h
    T = Ts*T0/369.89
    return (T,vl,vv)
end

#Optimization: does not calculate crit_pure for each value.
function wilson_k_values(model::MultiFluidModel,p,T,crit = nothing)
    n = length(model)
    K0 = zeros(typeof(p+T),n)
    pure = split_model.(model)
    _Tc = model.properties.Tc.values
    _Pc = model.properties.pc.values
    for i ∈ 1:n
        pure_i = pure[i]
        Tc,pc = _Pc[i],_Tc[i]
        ps = first(saturation_pressure(pure_i,0.7*Tc))
        ω = -log10(ps/pc) - 1.0
        K0[i] = exp(log(pc/p)+5.373*(1+ω)*(1-Tc/T))
    end
    return K0
end