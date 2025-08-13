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
    f(∂p) = simple_volume(model,∂p,T,z)
    V,∂V∂p = Solvers.f∂f(f,p)
    return SVector(V,∂V∂p)
end

function V∂V∂T(model,p,T,z::AbstractVector=SA[1.0])
    f(∂T) = simple_volume(model,p,∂T,z)
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
    ∂G∂T(x) = Solvers.derivative(G,x)
    ∂²G∂T²(x) = Solvers.derivative(∂G∂T,x)
    return ∂²G∂T²(T)
end
#property logic

function PT_property(model::GibbsBasedModel,p,T,z,phase,threaded,vol0,f::F,USEP::Val{UseP}) where {F,UseP}
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
    return -dvdT/v
end

function PT_property_gibbs(model,p,T,z,f::typeof(VT_isentropic_compressibility))
    ∂²g,∂g,g = ∂2𝕘(model,p,T,z)
    ∂²g∂T² = ∂²g[2,2]
    ∂²g∂p² = ∂²g[1,1]
    ∂²g∂T∂p = ∂²g[1,2]
    V = ∂g[1]
    return (∂²g∂T∂p*∂²g∂T∂p - ∂²g∂T²*∂²g∂p²)/∂²g∂T²/V
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
