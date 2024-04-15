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

    function M.component_fugacity_coefficient(model::C.EoSModel, cond, i, Z = C.compressibility_factor(model,cond.p,cond.T,cond.z), forces = nothing, s_v = nothing)
        p,T,z = cond.p,cond.T,cond.z
        R = C.Rgas(model)
        RT = R*T
        function fun(x)
            V = Z*C.Rgas(model)*T/p
            C.eos_res(model,V,T,x)
        end
        TT = eltype(p+T+first(z) + Z + one(eltype(model)))
        μᵢ = C.Solvers.grad_at_i(fun,z,i,TT)
        lnϕᵢ = μᵢ/RT - log(Z)
        return lnϕᵢ
    end

    function M.component_fugacity(model::C.EoSModel, cond, i, Z = C.compressibility_factor(model,cond.p,cond.T,cond.z), forces = nothing, s_v = nothing)
        lnϕᵢ = component_fugacity_coefficient(eos, cond, i, Z, forces, scalars)
        return exp(lnϕᵢ)*p*z[i]
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


    include("MultiComponentFlash/stability.jl")
    include("MultiComponentFlash/flash.jl")
    include("MultiComponentFlash/flow_coupler.jl")
    include("MultiComponentFlash/viscosity.jl")
    include("MultiComponentFlash/Clapeyron_flash.jl")
end #module
