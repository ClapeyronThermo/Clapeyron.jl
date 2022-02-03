struct LJRefParam <: EoSParam
    epsilon::PairParam{Float64}
    sigma::PairParam{Float64}
    Mw::SingleParam{Float64}
end

struct LJRefConsts <: EoSParam
    n::Vector{Float64}
    t::Vector{Float64}
    d::Vector{Int}
    c::Vector{Int}
    beta::Vector{Float64}
    gamma::Vector{Float64}
    eta::Vector{Float64}
    epsilon::Vector{Float64}
    function LJRefConsts()
        n = [
        0.52080730e-2,
        0.21862520e+1,
        -0.21610160e+1,
        0.14527000e+1,
        -0.20417920e+1,
        0.18695286e0,
        -0.90988445e-1,
        -0.49745610e0,
        0.10901431e0,
        -0.80055922e0,
        -0.56883900e0,
        -0.62086250e0,
        -0.14667177e+1,
        0.18914690e+1,
        -0.13837010e0,
        -0.38696450e0,
        0.12657020e0,
        0.60578100e0,
        0.11791890e+1,
        -0.47732679e0,
        -0.99218575e+1,
        -0.57479320e0,
        0.37729230e-2,
        ]
        
        t = [
        1.000
        0.320
        0.505
        0.672
        0.843
        0.898
        1.294
        2.590
        1.786
        2.770
        1.786
        1.205
        2.830
        2.548
        4.650
        1.385
        1.460
        1.351
        0.660
        1.496
        1.830
        1.616
        4.970
        ]
        
        d = [
        4
        1
        1
        2
        2
        3
        5
        2
        2
        3
        1
        1
        1
        1
        2
        3
        3
        2
        1
        2
        3
        1
        1
        ]
        c = 
        [
        1
        2
        1
        2
        2
        1
        
        ]
        
        g = [
        2.067 0.625 0.710 0.2053
        1.522 0.638 0.860 0.4090
        8.820 3.910 1.940 0.6000
        1.722 0.156 1.480 1.2030
        0.679 0.157 1.490 1.8290
        1.883 0.153 1.945 1.3970
        3.925 1.160 3.020 1.3900
        2.461 1.730 1.110 0.5390
        28.20 383.0 1.170 0.9340
        0.753 0.112 1.330 2.3690
        0.820 0.119 0.240 2.4300]
        
        eta = g[:,1]
        beta = g[:,2]
        gamma = g[:,3]
        epsilon = g[:,4]  
        return new(n,t,d,c,beta,gamma,eta,epsilon)
    end
    
end

is_splittable(::LJRefConsts) = false

struct LJRef <: EmpiricHelmholtzModel
    components::Vector{String}
    icomponents::UnitRange{Int}
    params::LJRefParam
    consts::LJRefConsts
    references::Vector{String}
end
@registermodel LJRef

idealmodel(x::LJRef) = x
export LJRef

function LJRef(components;
    userlocations=String[],
    verbose=false)
    params,sites = getparams(components, ["SAFT/PCSAFT"]; userlocations=userlocations, verbose=verbose)
    Mw = params["Mw"]
    params["sigma"].values .*= 1E-10
    sigma = sigma_LorentzBerthelot(params["sigma"])
    epsilon = epsilon_LorentzBerthelot(params["epsilon"], params["k"])
    params = LJRefParam(epsilon,sigma,Mw)
    consts = LJRefConsts()
    icomponents = 1:length(components)
    references = ["10.1063/1.4945000"]
    return LJRef(components,icomponents,params,consts,references)
end

function _f0(model::LJRef,ρ,T,z=SA[1.0],∑z = sum(z))
    ϵ = model.params.epsilon.values
    σ = model.params.sigma.values
    lnΣz = log(∑z)
    res = zero(ρ+T+first(z))
    for i  ∈ @comps
        τᵢ = 1.32/(T/ϵ[i])  
        δᵢ = (N_A*ρ*σ[i]^3)/0.31
        aᵢ = log(δᵢ) + 1.5*log(τᵢ) + 1.515151515*τᵢ + 6.262265814 
        res += z[i]*(aᵢ + log(z[i]) - lnΣz)
    end
    return res
end

function _fr(model::LJRef,δ,τ)
    ai = zero(δ+τ)
    n = model.consts.n
    t = model.consts.t
    d = model.consts.d
    c = model.consts.c
    β = model.consts.beta
    γ = model.consts.gamma
    η = model.consts.eta
    ε = model.consts.epsilon   

    @inbounds begin
        for k ∈ 1:6
            ai += n[k]*(δ^d[k])*(τ^t[k])
        end
        for (k,k_) ∈ zip(7:12,1:6)
            ai += n[k]*(δ^d[k])*(τ^t[k])*exp(-δ^c[k_])
        end

        for (k,k_) ∈ zip(13:23,1:11)
            ai += n[k]*(δ^(d[k]))*(τ^(t[k]))*
            exp(-η[k_]*(δ - ε[k_])^2 - β[k_]*(τ -γ[k_])^2)
        end
    end
    return ai
end

#TODO: better relations? EoSRef was done with one fluid in mind.
#this is technically an unsafe extension.
function _v_scale(model::LJRef,z=SA[1.0])
    σᵢᵢ = model.params.sigma.diagvalues
    val = N_A*sum(z[i]*σᵢᵢ[i]^3 for i in 1:length(z))
    return val
end

function T_scale(model::LJRef,z=SA[1.0])
    ϵ = model.params.epsilon.values
    σ = model.params.sigma.values
    σϵ_mix = zero(eltype(z))
    σ_mix = zero(eltype(z))
    comps = length(model)
    for i ∈ 1:comps
        zi = z[i]
        zii = zi*zi
        σ3 = σ[i,i]^3
        σϵ_mix += zii*σ3*ϵ[i,i]
        σ_mix += zii*σ3
        for j ∈ 1:(i-1)
            σ3ij = σ[i,j]^3
            zij = zi*z[j]
            σϵ_mix += 2*zij*σ3ij*ϵ[i,j]
            σ_mix += 2*zij*σ3ij*ϵ[i,j]
        end
    end
    return σϵ_mix/σ_mix
end

function lb_volume(model::LJRef, z = SA[1.0])
    val = π/6*_v_scale(model,z)
    return val
end

function eos(model::LJRef,V,T,z = SA[1.0])
    Σz = sum(z)
    ρ = Σz/V
    α0 = _f0(model,ρ,T,Σz)
    τ = 1.32/(T/T_scale(model,z))
    δ = (ρ*_v_scale(model,z)/Σz)/0.31
    αr =  _fr(model,δ,τ)
    x1 = R̄*T*Σz*αr 
    x2 =  R̄*T*α0
    return x1+x2
end

function a_ideal(model::LJRef,V,T,z = SA[1.0])
    Σz = sum(z)
    ρ = Σz/V
    α0 = _f0(model,ρ,T,Σz)
    return α0/Σz
end

function a_res(model::LJRef,V,T,z = SA[1.0])
    Σz = sum(z)
    ρ = Σz/V
    α0 = _f0(model,ρ,T,Σz)
    τ = 1.32/(T/T_scale(model,z))
    δ = (ρ*_v_scale(model,z)/Σz)/0.31
    return  _fr(model,δ,τ)
end

function eos_res(model::LJRef,V,T,z = SA[1.0])
    Σz = sum(z)
    ρ = Σz/V
    α0 = _f0(model,ρ,T,Σz)
    τ = 1.32/(T/T_scale(model,z))
    δ = (ρ*_v_scale(model,z)/Σz)/0.31
    αr =  _fr(model,δ,τ)
    return R̄*T*Σz*αr 
end
#=
function ljref_psat(Tr,pc)
    n = (0.54000e+1,0.44704e01,-0.18530e+1,0.19890e0,-0.11250e+1)
    t = (1.,1.5,4.7,2.5,21.4)  
    tr1 = one(Tr) - Tr
    res = sum(ni*tr1^ti for (ni,ti) in zip(n,t))
    return exp(res/Tr)*pc 
end
=#
function ljref_rholsat(Tr)
    n = (0.1362e+1,0.2093e+1,-0.2110e+1,0.3290e0,0.1410e+1)
    t = (0.313 ,0.940,1.630,17.,2.4)
    tr1 = one(Tr) - Tr
    res = sum(ni*tr1^ti for (ni,ti) in zip(n,t))
    rhoc = 0.31
    return (1+res)*rhoc
end

function ljref_rhovsat(Tr)
    n = (-0.69655e+1,-0.10331e+3,-0.20325e+1,-0.44481e+2,-0.18463e+2,-0.26070e+3)
    t = (1.320 ,19.24,0.360,8.780,4.040,41.60)
    tr1 = one(Tr) - Tr
    res = sum(ni*tr1^ti for (ni,ti) in zip(n,t))
    rhoc = 0.31
    return exp(res)*rhoc
end

function x0_sat_pure(model::LJRef,T)
    σ = model.params.sigma.values[1]
    Tc = 1.32*T_scale(model)
    ρl =  ljref_rholsat(T/Tc)/(N_A*σ^3)
    ρv =  ljref_rhovsat(T/Tc)/(N_A*σ^3)
    return [log10(1/ρl),log10(1/ρv)]
end

function p_scale(model::LJRef,z = SA[1.0])
    rhoc = 1/(_v_scale(model,z))
    Tc = T_scale(model,z)
    return R̄*Tc*rhoc
end



 