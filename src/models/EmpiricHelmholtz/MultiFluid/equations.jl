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
        mix_geomean,
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
    vc = model.properties.Vc.values
    #isone(length(z)) && return only(vc) 
    res = mixing_rule_asymetric(
        mix_mean3,
        _gerg_asymetric_mix_rule,
        z,
        vc,
        model.ideal.gamma_v.values,
        model.ideal.beta_v.values,
    )
    return res/(Σz*Σz)
end

function lb_volume(model::MultiFluidModel,z=SA[1.])
    return dot(z,model.properties.lb_volume.values)
end

function _delta(model::MultiFluidModel, rho, T, z=SA[1.],Σz = sum(z))
    vcmix = _v_scale(model,z,Σz)
    return rho * vcmix
end

function _tau(model::MultiFluidModel, rho, T, z=SA[1.],Σz = sum(z))
    Tcmix  = _T_scale(model,z,Σz)
    return Tcmix / T
end

function reduced_a_ideal(model::MultiFluidModel, ρ, T, z=SA[1.], Σz = sum(z))
    RR = 8.314472 / 8.314510
    #common_type = promote_type(typeof(ρ), typeof(T), eltype(x))
    _0 = zero(ρ + T + first(z))
    lnΣz = log(Σz)
    res = _0
    vc = model.properties.Vc.values
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

function reduced_a_res(model::MultiFluidModel,δ,τ,z)
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
        n_pol = view(nᵢ,k1)
        t_pol = view(tᵢ,k1)
        d_pol = view(dᵢ,k1)
        ai += term_ar_pol(δ,τ,lnδ,lnτ,_0,n_pol,t_pol,d_pol)
        
        n_exp = view(nᵢ,k2)
        t_exp = view(tᵢ,k2)
        d_exp = view(dᵢ,k2)
        c_exp = view(cᵢ,kexp)
        ai += term_ar_exp(δ,τ,lnδ,lnτ,_0,n_exp,t_exp,d_exp,c_exp,FillArrays.Ones(k2))
 
        res += z[i]*ai 
    end
    return res
end

function reduced_a_departure(model::MultiFluidModel,δ,τ,z)
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
            #GERG2008 uses a very particular set of terms
            #instead of -η(δ-ε)^2 - β(τ-γ)^2
            #uses -η(δ-ε)^2 - β(δ-γ)
            k1,k2,kgerg = ith_index(k_all,k_exp,ii)
            
            n_pol = view(nᵢⱼ,k1)
            t_pol = view(tᵢⱼ,k1)
            d_pol = view(dᵢⱼ,k1)
            aij += term_ar_pol(δ,τ,lnδ,lnτ,_0,n_pol,t_pol,d_pol)
            
            n_gauss = view(nᵢⱼ,k2)
            t_gauss = view(tᵢⱼ,k2)
            d_gauss = view(dᵢⱼ,k2)
            η = view(ηᵢⱼ,kgerg)
            β = view(βᵢⱼ,kgerg)
            γ = view(γᵢⱼ,kgerg)
            ε = view(εᵢⱼ,kgerg)
            aij += term_ar_gerg2008(δ,τ,lnδ,lnτ,_0,n_gauss,t_gauss,d_gauss,η,β,γ,ε)
            
           res +=z[i]*z[j]*Fᵢⱼ*aij
        end
     end
    return res
end

function a_ideal(model::MultiFluidModel, V, T, z=SA[1.0])
    Σz = sum(z)
    ρ = Σz*1.0e-3/V
    reduced_a_ideal(model, ρ, T, z,Σz)/Σz
end

function a_res(model::MultiFluidModel, V, T, z=SA[1.0])
    Σz = sum(z)
    invn = 1/Σz
    invn2 = invn*invn
    ρ = Σz*1.0e-3/V
    δ = _delta(model, ρ, T, z,Σz)
    τ = _tau(model, ρ, T, z,Σz)
    return (reduced_a_res(model,δ,τ,z)*invn + reduced_a_departure(model,δ,τ,z)*invn2)
end

function eos(model::MultiFluidModel, V, T, z=SA[1.0])
    Σz = sum(z)
    invn = 1/Σz
    ρ = Σz*1.0e-3/V
    δ = _delta(model, ρ, T, z,Σz)
    τ = _tau(model, ρ, T, z,Σz)
    return  R̄*T*(reduced_a_ideal(model,ρ,T,z,Σz) + reduced_a_res(model,δ,τ,z)+reduced_a_departure(model,δ,τ,z)*invn)
end

function eos_res(model::MultiFluidModel, V, T, z=SA[1.0])
    Σz = sum(z)
    invn = 1/Σz
    ρ = Σz*1.0e-3/V
    δ = _delta(model, ρ, T, z,Σz)
    τ = _tau(model, ρ, T, z,Σz)
    return  R̄*T*(reduced_a_res(model,δ,τ,z)+reduced_a_departure(model,δ,τ,z)*invn)
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
        Tc,pc = _Tc[i],_Pc[i]
        ps = first(saturation_pressure(pure_i,0.7*Tc))
        ω = -log10(ps/pc) - 1.0
        K0[i] = exp(log(pc/p)+5.373*(1+ω)*(1-Tc/T))
    end
    return K0
end
