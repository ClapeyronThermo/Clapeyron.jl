module ClapeyronMultiComponentFlashExt
    using Clapeyron: Clapeyron
    using MultiComponentFlash: MultiComponentFlash
    const C = Clapeyron
    const M = MultiComponentFlash
    const S = C.StaticArrays
    const ForwardDiff = C.ForwardDiff



    ##
    ##passing GenericCubicEoS to Clapeyron.jl userlocations
    ##

    C.can_nt(x::M.GenericCubicEOS) = true
    C.can_nt(x::M.MultiComponentMixture) = true

    function C.to_nt(x::M.MultiComponentMixture)
        res = Dict{Symbol,Any}()
        if x.binary_interaction !== nothing
            res[:k] = x.binary_interaction
        end
        props = x.properties
        res[:Mw] = [1000*getfield(props_i,:mw) for props_i in props]
        res[:Tc] = [getfield(props_i,:T_c) for props_i in props]
        res[:Pc] = [getfield(props_i,:p_c) for props_i in props]
        res[:Vc] = [getfield(props_i,:V_c) for props_i in props]
        res[:acetricfactor] = [getfield(props_i,:ω) for props_i in props]
        return res
    end
    function C.to_nt(x::M.GenericCubicEOS)
        res = C.to_nt(x.mixture)
        #TODO:check that volume shift
    end

    ##
    ##EoS properties to MultiComponentFlash.jl properties
    ##
    M.number_of_components(model::C.EoSModel) = length(model)

    function M.mixture_compressibility_factor(eos::C.EoSModel, cond,
        forces = M.force_coefficients(eos, cond),
        scalars = M.force_scalars(eos, cond, forces))
        phase::Symbol = get(cond,:phase,:unknown)
        C.compressibility_factor(eos,cond.p,cond.T,cond.z,phase = phase)
    end

    function M.initial_guess_K!(K, eos::C.EoSModel, cond)
        C.tp_flash_K0!(K, eos, cond.p, cond.T)
    end

    #this is only defined with cubic EoS.
    M.force_coefficients(eos::C.EoSModel, cond;static_size = false) = nothing
    M.force_scalars(eos::C.EoSModel, cond, forces) = nothing
    M.force_coefficients!(forces, eos::C.EoSModel, c) = nothing

    function M.mixture_fugacities!(f, eos::C.EoSModel, cond, forces = M.force_coefficients(eos, cond), scalars = M.force_scalars(eos, cond, forces))
        phase::Symbol = get(cond,:phase,:unknown)
        C.fugacity_coefficient!(f,eos,cond.p,cond.T,cond.z,phase = phase)
        f .*= cond.p
        f .= f .* cond.z
        return f
    end

    function M.component_fugacity(model::C.EoSModel, cond, i, Z = C.compressibility_factor(model,cond.p,cond.T,cond.z), forces = nothing, s_v = nothing)
        p,T,z = cond.p,cond.T,cond.z
        R = C.Rgas(model)
        RT = R*T
        function fun(x)
            V = Z*C.Rgas(model)*T/p
            C.eos_res(model,V,T,x)
        end
        TT = eltype(p+T+first(z) + Z + one(eltype(model)))
        μᵢ = C.Solvers.grad_at_i(fun,z,i,TT)
        ϕᵢ = exp(μᵢ/RT)/Z
        return ϕᵢ*p*z[i]
    end
    if isdefined(M,:eostype)
        M.eostype(model::C.EoSModel) = Base.summary(model)
    end

    if isdefined(M,:molar_masses)
        M.molar_masses(model::C.EoSModel) = C.mw(model) .* 0.001
    end

    if isdefined(M,:component_names)
        M.component_names(model::C.EoSModel) = model.components #TODO: we need an abstraction of this type on our code
    end

    function M.lbc_viscosity(eos::C.EoSModel, p, temperature, ph::M.FlashedPhase{T}; coeff = (0.1023, 0.023364, 0.058533, -0.040758, 0.0093324), shift = -1e-4) where T
        z = ph.mole_fractions
        
        mw_mix, P_pc, T_pc, V_pc = M.pseudo_critical_properties(eos, z)
        mu_atm = M.atmospheric_mu_estimate(eos, z, temperature)
        e_mix = M.mixture_viscosity_parameter(mw_mix, T_pc, P_pc)
        # From Jossi et al
        # coeffs = [0.1023, 0.023364, 0.058533, -0.040758, 0.0093724, 1e-4]
        # From LBC paper
        # coeffs = [0.1023, 0.023364, 0.058533, -0.040758, 0.0093324, -1e-4]
        # Compute reduced density
        V = M.molar_volume(eos, p, temperature, ph)
        rho_r = V_pc/V
        corr = zero(T)
        for (i, c) in enumerate(coeff)
            corr += c*rho_r^(i-1)
        end
        mu_correction = (corr^4 + shift)/e_mix
        # Final expression is compound and given in centi poise. We convert to Pa s
        # instead for strict SI outputs.
        mu = 1e-3*(mu_atm + mu_correction)
        return mu::T
    end

    function M.pseudo_critical_properties(model::C.EoSModel, z::V) where V<:AbstractVector{T} where T
        mw = C.molecular_weight(model,z)
        Tc = C.dot(model.params.Tc.values,z)
        Pc = C.dot(model.params.Pc.values,z)
        Vc = C.volume(model,Pc,Tc,z) #we calculate the critical volume of the EoS instead of doing it component wise?
        return mw,Pc,Tc,Vc 
    end

    function M.atmospheric_mu_estimate(model::C.EoSModel, z::V, temperature) where V<:AbstractVector{T} where T
        a = zero(T)
        b = zero(T)
        Tc = model.params.Tc.values
        Pc = model.params.Pc.values
        Mw = C.mw(model)
        @inbounds for i in 1:length(model)
            mw = Mw[i]
            T_c = Tc[i]
            p_c = Pc[i]
            zi = z[i]
            # Add contributions to atmospheric mu
            T_r = temperature/T_c
            e_i = M.mixture_viscosity_parameter(mw, T_c, p_c)
            if T_r > 1.5
                mu_i = 17.78e-5*(4.58*T_r - 1.67)^0.625
            else
                mu_i = 34e-5*T_r^(0.94)
            end
            tmp = sqrt(1000.0*mw)*zi
            a += tmp*mu_i/e_i
            b += tmp
        end
        return a/b
    end
    
    include("MultiComponentFlash/stability.jl")
    include("MultiComponentFlash/flash.jl")
    include("MultiComponentFlash/flow_coupler.jl")
    include("MultiComponentFlash/Clapeyron_flash.jl")
end #module
