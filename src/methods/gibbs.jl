#abstract type GibbsBasedModel <: EoSModel end

#derivative logic

function ∂𝕘∂T(model,p,T,z::AbstractVector)
    g(∂T) = eos_g(model,p,∂T,z)
    return Solvers.derivative(g,T)
end

function ∂𝕘∂p(model,p,T,z::AbstractVector)
    g(∂p) = eos_g(model,∂p,T,z)
    return Solvers.derivative(g,p)
end

function ∂𝕘(model,p,T,z)
    f(∂p,∂T) = eos_g(model,∂p,∂T,z)
    _f,_df = Solvers.fgradf2(f,p,T)
    return _df,_f
end

function ∂𝕘_vec(model,p,T,z::AbstractVector)
    _df,_f = ∂𝕘(model,p,T,z)
    return SVector(_f,_df[1],_df[2])
end

function 𝕘∂𝕘dp(model,p,T,z::AbstractVector)
    f(x) = eos_g(model,x,T,z)
    G,∂G∂p = Solvers.f∂f(f,p)
    return SVector(G,∂G∂p)
end

function 𝕘∂𝕘dT(model,p,T,z::AbstractVector)
    f(x) = eos_g(model,p,x,z)
    G,∂G∂T = Solvers.f∂f(f,T)
    return SVector(G,∂G∂T)
end

function V∂V∂p(model,p,T,z::AbstractVector=SA[1.0])
    f(∂p) = ∂𝕘∂p(model,∂p,T,z)
    V,∂V∂p = Solvers.f∂f(f,p)
    return SVector(V,∂V∂p)
end

function V∂V∂T(model,p,T,z::AbstractVector=SA[1.0])
    f(∂T) = ∂𝕘∂p(model,p,∂T,z)
    V,∂V∂T = Solvers.f∂f(f,T)
    return SVector(V,∂V∂T)
end

function ∂2𝕘(model,p,T,z)
    f(_p,_T) = eos_g(model,_p,_T,z)
    _f,_∂f,_∂2f = Solvers.∂2(f,p,T)
    return (_∂2f,_∂f,_f)
end

function 𝕘_hess(model,p,T,z)
    f(w) = eos_g(model,first(w),last(w),z)
    p,T = promote(p,T)
    pT_vec = SVector(p,T)
    return Solvers.hessian(f,pT_vec)
end

function ∂²𝕘∂T²(model,p,T,z)
    G(x) = eos_g(model,p,x,z)
    _,_,∂²G∂T² = Solvers.f∂f∂2f(G,T)
    return ∂²G∂T²
end
#property logic

function PT_property(model::GibbsBasedModel,p,T,z,phase,threaded,vol0,f::F,vol) where {F}
    z isa Number && return PT_property_gibbs(model,p,T,SVector(z),f)
    return PT_property_gibbs(model,p,T,z,f)
end

PT_property_gibbs(model,p,T,z,f::typeof(pressure)) = p
PT_property_gibbs(model,p,T,z,f::typeof(temperature)) = T
PT_property_gibbs(model,p,T,z,f::typeof(volume)) = volume(model,p,T,z)
PT_property_gibbs(model,p,T,z,f::typeof(VT_entropy)) = -∂𝕘∂T(model,p,T,z)
PT_property_gibbs(model,p,T,z,f::typeof(VT_gibbs_free_energy)) = eos_g(model,p,T,z)
PT_property_gibbs(model,p,T,z,f::typeof(VT_helmholtz_free_energy)) = eos_g(model,p,T,z) - p*volume(model,p,T,z)

PT_property_gibbs(model,p,T,z,f::typeof(VT_mass_entropy)) = PT_property_gibbs(model,p,T,z,VT_entropy)/molecular_weight(model,z)
PT_property_gibbs(model,p,T,z,f::typeof(VT_mass_gibbs_free_energy)) = PT_property_gibbs(model,p,T,z,VT_gibbs_free_energy)/molecular_weight(model,z)
PT_property_gibbs(model,p,T,z,f::typeof(VT_mass_helmholtz_free_energy)) = PT_property_gibbs(model,p,T,z,VT_helmholtz_free_energy)/molecular_weight(model,z)
PT_property_gibbs(model,p,T,z,f::typeof(VT_mass_isobaric_heat_capacity)) = PT_property_gibbs(model,p,T,z,VT_isobaric_heat_capacity)/molecular_weight(model,z)
PT_property_gibbs(model,p,T,z,f::typeof(VT_mass_isochoric_heat_capacity)) = PT_property_gibbs(model,p,T,z,VT_isochoric_heat_capacity)/molecular_weight(model,z)

volume_impl(model::GibbsBasedModel,p,T,z,phase,threaded,vol0) = ∂𝕘∂p(model,p,T,z)

function PT_property_gibbs(model,p,T,z,f::typeof(VT_internal_energy))
    g,∂g∂p,∂g∂T = ∂𝕘_vec(model,p,T,z)
    return g - T*∂g∂T - p*∂g∂p
end

function PT_property_gibbs(model,p,T,z,f::typeof(VT_enthalpy))
    g,∂g∂T = 𝕘∂𝕘dT(model,p,T,z)
    return g - T*∂g∂T
end

function PT_property_gibbs(model,p,T,z,f::typeof(VT_isobaric_heat_capacity))
    ∂²g∂T² = ∂²𝕘∂T²(model,p,T,z)
    return -T*∂²g∂T²
end

function PT_property_gibbs(model,p,T,z,f::typeof(VT_isochoric_heat_capacity))
    ∂²g = 𝕘_hess(model,p,T,z)
    ∂²g∂T² = ∂²g[2,2]
    ∂²g∂p² = ∂²g[1,1]
    ∂²g∂T∂p = ∂²g[1,2]
    return -T*(∂²g∂T² + ∂²g∂p²/∂²g∂T∂p)
end

function PT_property_gibbs(model,p,T,z,f::typeof(VT_adiabatic_index))
    ∂²g = 𝕘_hess(model,p,T,z)
    ∂²g∂T² = ∂²g[2,2]
    ∂²g∂p² = ∂²g[1,1]
    ∂²g∂T∂p = ∂²g[1,2]
    Cv = -T*(∂²g∂T² + ∂²g∂p²/∂²g∂T∂p)
    Cp = -T*∂²g∂T²
    return 1/(1 + ∂²g∂p²/(∂²g∂T∂p*∂²g∂T²))
end

function PT_property_gibbs(model,p,T,z,f::typeof(VT_isothermal_compressibility))
    V,∂V∂p = V∂V∂p(model,p,T,z)
    return -∂V∂p/V
end

function PT_property_gibbs(model,p,T,z,f::typeof(VT_isobaric_expansivity))
    v,dvdT = V∂V∂T(model,p,T,z)
    return dvdT/v
end

function PT_property_gibbs(model,p,T,z,f::typeof(VT_isentropic_compressibility))
    ∂²g,∂g,g = ∂2𝕘(model,p,T,z)
    ∂²g∂T² = ∂²g[2,2]
    ∂²g∂p² = ∂²g[1,1]
    ∂²g∂T∂p = ∂²g[1,2]
    V = ∂g[1]
    return (∂²g∂T∂p*∂²g∂T∂p - ∂²g∂T²*∂²g∂p²)/∂²g∂T²/V
end

function PT_property_gibbs(model,p,T,z,f::typeof(VT_speed_of_sound))
    Mr = molecular_weight(model,z)
    ∂²g,∂g,g = ∂2𝕘(model,p,T,z)
    ∂²g∂T² = ∂²g[2,2]
    ∂²g∂p² = ∂²g[1,1]
    ∂²g∂T∂p = ∂²g[1,2]
    V = ∂g[1]
    βsρ = (∂²g∂T∂p*∂²g∂T∂p - ∂²g∂T²*∂²g∂p²)/∂²g∂T²
    V*sqrt(1/(βsρ*Mr))
end

function VT_pressure(model::GibbsBasedModel,V,T,z)
    _p0 = x0_pressure(model,V,T,z)
    p0 = one(Base.promote_eltype(model,V,T,z))*_p0
    !isfinite(p0) && return p0
    function fixpoint_p(pᵢ)
        Vᵢ,∂V∂pᵢ = V∂V∂p(model,pᵢ,T,z)
        dp = log(V/Vᵢ)*Vᵢ/∂V∂pᵢ
        px = pᵢ + dp
        return pᵢ + dp
    end
    return Solvers.fixpoint(fixpoint_p,p0,Solvers.SSFixPoint(),rtol = 1e-12)
end

x0_volume_liquid(model::GibbsBasedModel,p,T,z) = simple_volume(model,p,T,z)
x0_volume_solid(model::GibbsBasedModel,p,T,z) = simple_volume(model,p,T,z)
x0_volume_gas(model::GibbsBasedModel,p,T,z) = simple_volume(model,p,T,z)

function x0_pressure(model,V,T,z)
    p = p_scale(model,z)*one(T+first(z)+V)
        for i in 1:20
        if volume(model,p,T) <= V || !isfinite(p)
            return p
        end
        p *= 2
    end
    return p
end

function chemical_potential_impl(model::GibbsBasedModel,p,T,z,phase,threaded,vol0)
    return VT_molar_gradient(model,p,T,z,eos_g)
end

function __calculate_reference_state_consts(model::GibbsBasedModel,v,T,p,z,H0,S0,phase)
    S00 = entropy(model,p,T,z)
    a1 = (S00 - S0)#/∑z
    H00 = enthalpy(model,p,T,z)
    a0 = (-H00 + H0)#/∑z
    return a0,a1
end

"""
    type,p,T,W = gibbsmodel_reference_state_consts(model)
    type,p,T,W = gibbsmodel_reference_state_consts(model,other_model)
    
Returns an equilibrium condition to equilibrate the Gibbs energies of two models.
Used for solid-fluid equilibria.
By default, it returns `nothing`. 
The two-argument method is used to disambiguate between two different models.
Available options for the type are:
    - :dH: difference in enthalpy at p,T conditions (`h(other_model) - h(model)`) is equal to W
    - :zero: the models are already equilibrated, no additional calculation is necessary (like `IAPWS06` in conjunction with `IAPWS05`)


The equilibration corresponds to the calculation of constants `k1` and `k2`, that enforce the Gibbs criteria: `gibbs_energy(model,p,T) + k1 + k2*T == gibbs_energy(other_model,p,T)`
The constants `k1` and `k2` are calculated by `Clapeyron.calculate_gibbs_reference_state(model,other_model)`
"""
gibbsmodel_reference_state_consts(model::EoSModel) = nothing
gibbsmodel_reference_state_consts(model1,model2) = nothing
function _gibbsmodel_reference_state_consts(model1,model2)
    ref1 = gibbsmodel_reference_state_consts(model1,model2)
    ref2 = gibbsmodel_reference_state_consts(model2,model1)
    ref2 == nothing && ref1 == nothing && return nothing,0
    ref1 == nothing && return ref2,2
    ref2 == nothing && return ref1,1
    throw(error("invalid specification for gibbs_reference_state_consts: both model1 and model2 define their own order."))
end


"""
    k1,k2 = calculate_gibbs_reference_state(model,other_model)
    
Calculates the reference state constants that forces the equilibrium conditions specified by `Clapeyron.gibbsmodel_reference_state_consts`.
"""
function calculate_gibbs_reference_state(model1::EoSModel,model2::EoSModel,x1 = SA[1.0],x2 = SA[1.0])

    _0 = zero(Base.promote_eltype(model1,model2))
    (model1 isa GibbsBasedModel) || (model2 isa GibbsBasedModel) || return _0,_0 
    refx,nx = _gibbsmodel_reference_state_consts(model1,model2)
    if refx == nothing
        ref1 = gibbsmodel_reference_state_consts(model1)
        ref2 = gibbsmodel_reference_state_consts(model2)
        if ref1 == nothing && ref2 == nothing
            throw(error("Empty Gibbs reference. For Gibbs models, define `Clapeyron.gibbsmodel_reference_state_consts(model)`"))
        end
        if ref1 == nothing
            ref = ref2
            n = 2
        elseif ref2 == nothing
            ref = ref1
            n = 1
        elseif ref2 != nothing && ref1 != nothing
            
            isnothing(refx) && throw(error("Empty Gibbs reference. For Gibbs models, define `Clapeyron.gibbsmodel_reference_state_consts(model1,model2)`"))
        end
    else
        ref = refx
        n = nx
    end
    type,p,T,W = ref

    if type == :dH
        #=
        dH reference: at p = p0,T = T0, g1 = g2, H1 - H2 = Hfus:
        eos_g(model1) + g1 + g2*T0 = eos_g(model2)
        enthalpy(model2) - enthalpy(model1) - W - g1 = 0
        =#
        gibbs1,gibb2 = gibbs_energy(model1,p,T,x1),gibbs_energy(model2,p,T,x2)
        h1,h2 = enthalpy(model1,p,T,x1),enthalpy(model2,p,T,x2)
        dH = W
        if n == 2 #invert
            gibbs1,gibb2 = gibbs2,gibbs1
            h1,h2 = h2,h1
            dH = -dH
        end
        g1 = h1 - h2 + dH
        g2 = (gibb2 - g1 - gibbs1)/T
        return g1,g2
    elseif type == :zero
        return _0,_0
    else
        throw(error("invalid Gibbs reference state. Expected :dH, got: $type"))
    end
end



function dpdT_saturation_gibbs(model1,model2,p,T,w1 = SA[1.0],w2 = SA[1.0],k = calculate_gibbs_reference_state(model1,model2);phase1 = :unknown,phase2 = :unknown)
    k1,k2 = k
    ∑w1 = sum(w1)
    ∑w2 = sum(w2)  
    v1 = volume(model1,p,T)/∑w1
    v2 = volume(model2,p,T)/∑w2
    s1 = entropy(model1,p,T,w1)
    s2 = entropy(model2,p,T,w2)
    dS = s1/∑w1 + k2 - s2/∑w2
    dv = v1/∑w1 - v2/∑w2
    return dS/dv
end

#=
init of pressure-based iterative methods
=#

function gibbs2_expansion(model::GibbsBasedModel,p,T)
    f(_p) = gibbs_energy(model,_p,T)
    return Solvers.f∂f∂2f(f,p)
end

function gibbs2_expansion(model,p,T)
    V = volume(model,p,T)
    f(_V) = eos(model,_V,T)
    a,da,d2a = Solvers.f∂f∂2f(f,V)
    g = a + p*V
    dg = V
    d2g = -1/d2a
    return g,dg,d2g
end

function g_and_v(model,p,T,v;phase = :unknown)
    v = volume(model,p,T,SA[1.0],phase = phase,vol0 = v)
    g = VT_gibbs_free_energy(model,v,T,SA[1.0])
    return g,v
end

function g_and_v(model::GibbsBasedModel,p,T,v;phase = :unknown)
    return 𝕘∂𝕘dp(model,p,T,SA[1.0])
end

function g_and_sv(model,p,T,v;phase = :unknown)
    v = volume(model,p,T,SA[1.0],phase = phase,vol0 = v)
    g = VT_gibbs_free_energy(model,v,T,SA[1.0])
    s = VT_entropy(model,v,T,SA[1.0])
    return g,s,v
end

function g_and_sv(model::GibbsBasedModel,p,T,v;phase = :unknown)
    g,v,sn =  ∂𝕘_vec(model,p,T,SA[1.0])
    return g,-sn,v
end

function eos_g_incomp(model,p,T,z,p0,T0)
    V = simple_volume(model,p,T,z)
    return V*(p - p0) + gibbs_cp_integral(idealmodel(model),T,z,T0)
end
