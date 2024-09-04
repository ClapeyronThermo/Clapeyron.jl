#this is a julia port of https://github.com/usnistgov/teqp/pull/154
#https://github.com/usnistgov/teqp/blob/main/include/teqp/algorithms/phase_equil.hpp

module PhaseEq
using ForwardDiff
using Clapeyron
using LinearAlgebra
using Clapeyron: viewn,viewfirst,viewlast
using Clapeyron.Solvers
using Clapeyron.DiffResults

mutable struct RequiredPhaseDerivatives{T}
    ρ::T
    R::T
    Ψᵣ::T
    ∇Ψᵣ::Vector{T}
    ∇²Ψᵣ::Matrix{T}
    ∂Ψᵣ∂T::T
    ∂₂Ψᵣ∂T²::T
    ∇∂Ψᵣ∂T::Vector{T}
    cache::DiffResults.MutableDiffResult{2, T, Tuple{Vector{T}, Matrix{T}}}
end

function RequiredPhaseDerivatives(model,T,ρᵢ,entropy = false)
    TT = Base.promote_eltype(model,T,ρᵢ)
    n = length(ρᵢ)
    result = RequiredPhaseDerivatives{TT}(n)
    update!(result,model,T,ρᵢ,entropy)
end

function RequiredPhaseDerivatives{T}(n) where T
    _0 = zero(T)
    _x = fill(_0,n)
    cache = DiffResults.HessianResult(_x)
    return RequiredCaloricDerivatives(_0,_0,_0,_x,fill(_0,(n,n)),_0,_0,fill(_0,n),cache)
end

function update!(m::RequiredPhaseDerivatives,model,T,ρᵢ)
    m.ρ = ρ = sum(ρᵢ)
    m.R = R = Rgas(model)
    ΨᵣT(_T) = Clapeyron.Ψ_res(model,T,ρᵢ)
    Ψᵣ,∂Ψᵣ∂T,∂₂Ψᵣ∂T² = Solvers.f∂f∂2f(_ΨᵣT,T)
    m.Ψᵣ = Ψᵣ
    m.∂Ψᵣ∂T = ∂Ψᵣ∂T
    m.∂₂Ψᵣ∂T² = ∂₂Ψᵣ∂T²
    ∇Ψ,HΨ = Clapeyron.Ψ_grad_and_hessian_res!(m.result,model,T,ρᵢ)
    m.∇Ψᵣ = ∇Ψ
    m.∇²Ψᵣ = HΨ
    Clapeyron.dΨdT_grad!(m.∇∂Ψᵣ∂T,model,T,ρᵢ)
    return result
end

#pressure
function p(res::RequiredPhaseDerivatives,T,ρᵢ)
    res.ρ*res.R*T - res.Ψᵣ + dot(res.∇Ψᵣ,ρᵢ)
end

function dpdT(m::RequiredPhaseDerivatives,T,ρᵢ)
    res.ρ*res.R - res.∂Ψᵣ∂T + dot(res.∇∂Ψᵣ∂T,ρᵢ)
end

function dpdρᵢ(m::RequiredPhaseDerivatives,T,ρᵢ)
    return dpdρᵢ!(zero(ρᵢ),m,T,ρᵢ)
end

function dpdρᵢ!(result,m::RequiredPhaseDerivatives,T,ρᵢ)
    mul!(result,transpose(ρᵢ),res.∇²Ψᵣ)
    result .+= res.R*T
    return result
end

#some additional notes
#apart from dpdρᵢ!, all similar mutating functions add the calculation to the result.
#this allows to avoid storing multiple intermediate results in memory
#setting the zeros is defered to the user, or you can use the nonmutating functions instead.

mutable struct RequiredCaloricDerivatives{T}
    ρ::T
    R::T
    Ψ₀::T
    ∂Ψ₀∂T::T
    ∂₂Ψ₀∂T²::T
    ∂₂Ψᵣ∂T²::T
    ∇Ψ₀::Vector{T}
    ∇∂Ψ₀∂T::Vector{T}
end

function RequiredCaloricDerivatives(model,T,ρᵢ,entropy = false)
    TT = Base.promote_eltype(model,T,ρᵢ)
    n = length(ρᵢ)
    result = RequiredCaloricDerivatives{TT}(n)
    update!(result,model,T,ρᵢ,entropy)
end

function RequiredCaloricDerivatives{T}(n = 0) where T
    _0 = zero(T)
    return RequiredCaloricDerivatives(_0,_0,_0,_0,_0,_0,fill(_0,n),fill(_0,n))
end

function update!(result::RequiredCaloricDerivatives,model::Clapeyron.EoSModel,T,ρᵢ,entropy = false)
    ∑z = sum(z)
    ρ = ∑z/V
    R = Clapeyron.Rgas(model)
    ig = Clapeyron.idealmodel(model)
    _Ψ₀T(_T) = Clapeyron.Ψ(ig,_T,ρᵢ)
    Ψ₀,∂Ψ₀∂T,∂₂Ψ₀∂T² = Solvers.f∂f∂2f(_Ψ₀T,T)
    if entropy
        _ΨᵣT(_T) = Clapeyron.Ψ_res(model,T,ρᵢ)
        _,_,∂₂Ψᵣ∂T² = Solvers.f∂f∂2f(_ΨᵣT,T)
    else
        ∂₂Ψᵣ∂T² = zero(Ψ₀)/zero(Ψ₀)
    end
    R = Rgas(model)
    result.R = R
    result.ρ = ρ
    result.Ψ₀ = Ψ₀
    result.∂Ψ₀∂T = ∂Ψ₀∂T
    result.∂₂Ψ₀∂T² = ∂₂Ψ₀∂T²
    result.∂₂Ψᵣ∂T² = ∂₂Ψᵣ∂T²
    Clapeyron.Ψ_grad!(result.∇Ψ₀,ig,T,ρᵢ)
    Clapeyron.dΨdT_grad!(result.∇∂Ψ₀∂T,model,T,ρᵢ)
    return result
end

#entropy
function s(ig::RequiredCaloricDerivatives,T,ρᵢ,res::RequiredPhaseDerivatives)
    return -1/ig.ρ * (ig.∂Ψ₀∂T + res.∂Ψᵣ∂T)
end

function dsdT(ig::RequiredCaloricDerivatives,T,ρᵢ,res::RequiredPhaseDerivatives)
    return -1/ig.ρ * (ig.∂₂Ψ₀∂T² + res.∂₂Ψᵣ∂T²)
end

function dsdρᵢ(ig::RequiredCaloricDerivatives,T,ρᵢ,res::RequiredPhaseDerivatives)
    return dsdρᵢ!(zero(ρᵢ),ig,T,ρᵢ,res)
end

function dsdρᵢ!(result,ig::RequiredCaloricDerivatives,T,ρᵢ,res::RequiredPhaseDerivatives)
    ρ = ig.ρ
    result .+= ig.∇∂Ψ₀∂T .+ res.∇∂ΨᵣT
    result .*= -1/ρ
    result .+= (1/ρ/ρ)*(ig.∂Ψ₀∂T + res.∂Ψᵣ∂T)
    return result
end

#helmholtz
function a(ig::RequiredCaloricDerivatives,T,ρᵢ,res::RequiredPhaseDerivatives)
    return (ig.Ψ₀ .+ res.Ψᵣ)/ig.ρ
end

function dadT(ig::RequiredCaloricDerivatives,T,ρᵢ,res::RequiredPhaseDerivatives)
    (ig.∂Ψ₀∂T .+ res.∂Ψᵣ∂T)/ig.ρ
end

function dadρᵢ(ig::RequiredCaloricDerivatives,T,ρᵢ,res::RequiredPhaseDerivatives)
    return dadρᵢ!(zero(ρᵢ),ig,T,ρᵢ,res)
end

function dadρᵢ!(result,ig::RequiredCaloricDerivatives,T,ρᵢ,res::RequiredPhaseDerivatives)
    ρ = ig.ρ
    result .+= ig.∇Ψ₀ .+ res.∇Ψᵣ
    result .*= 1/ρ
    result .-= (1/ρ/ρ)*(ig.Ψ₀ + res.Ψᵣ)
    return result
end

#internal energy
function u(ig::RequiredCaloricDerivatives,T,ρᵢ,res::RequiredPhaseDerivatives)
    return a(ig,T,ρᵢ,res) + s(ig,T,ρᵢ,res)
end

function dudT(ig::RequiredCaloricDerivatives,T,ρᵢ,res::RequiredPhaseDerivatives)
    return dadT(ig,T,ρᵢ,res) + s(ig,T,ρᵢ,res) + T*dsdT(ig,T,ρᵢ,res)
end

function dudρᵢ(ig::RequiredCaloricDerivatives,T,ρᵢ,res::RequiredPhaseDerivatives)
    return dudρᵢ!(zero(ρᵢ),ig,T,ρᵢ,res)
end

function dudρᵢ!(result,ig::RequiredCaloricDerivatives,T,ρᵢ,res::RequiredPhaseDerivatives)
    dsdρᵢ!(result,ig,T,ρᵢ,res)
    result .*= T
    return dadρᵢ!(result,ig,T,ρᵢ,res)
end

#enthalpy
function h(ig::RequiredCaloricDerivatives,T,ρᵢ,res::RequiredPhaseDerivatives)
    return u(ig,T,ρᵢ,res) + p(res,T,ρᵢ)(ig.ρ)
end

function dhdT(ig::RequiredCaloricDerivatives,T,ρᵢ,res::RequiredPhaseDerivatives)
    return dudT(ig,T,ρᵢ,res) + dpdT(res,T,ρᵢ)/ig.ρ
end

function dhdρᵢ(ig::RequiredCaloricDerivatives,T,ρᵢ,res::RequiredPhaseDerivatives)
    return dhdρᵢ!(zero(ρᵢ),ig,T,ρᵢ,res)
end

function dhdρᵢ!(result,ig::RequiredCaloricDerivatives,T,ρᵢ,res::RequiredPhaseDerivatives)
    ρ = ig.ρ
    dpdρᵢ!(result,res,T,ρᵢ)
    result .*= 1/ρ
    dudρᵢ(result,ig,T,ρᵢ,res)
    result .-= p(res,T,ρᵢ)/ρ/ρ
    return result
end

struct EqObject{M,T}
    model::M
    Nphases::Int
    Ncomponents::Int
    Nindependent::Int
    p0::T
    dpdT0::T
    dpdρᵢ0::Vector{T}
    derivatives::Vector{RequiredPhaseDerivatives{T}}
    caloric::Vector{RequiredCaloricDerivatives{T}}
    spec1::Specification{T}
    spec2::Specification{T}
end

_has_entropy(eq::EqObject) = spec1.sp == MolarEntropySpecification || spec2.sp == MolarEntropySpecification

function _has_caloric(eq::EqObject)
    sp1,sp2 = spec1.sp,spec2.sp
    caloric = (MolarEnthalpySpecification,MolarEntropySpecification,MolarInternalEnergySpecification)
    return (sp1 in caloric) | (sp2 in caloric)
end

function temperature(eq::EqObject,x)
    sp1,sp2 = spec1.sp,spec2.sp
    if sp1 == TemperatureSpecification
        return spec1.val
    elseif sp2 == TemperatureSpecification
        return spec2.val
    else
        return first(x)
    end
end

function update!(eq::EqObject,x)
    has_entropy = _has_entropy(eq)
    has_caloric = _has_caloric(eq)
    T = temperature(eq,x)
    x_rest = @view(x[2:end])
    derivatives = eq.derivatives
    caloric_derivatives = eq.caloric
    for i in 1:eq.Nphases
        ρᵢ = viewn(x_rest,eq.Ncomponents,i)
        update!(derivatives[i],model,T,ρᵢ)
        if has_caloric
            update!(caloric_derivatives[i],model,T,ρᵢ,has_entropy)
        end
    end
end

abstract type AbstractSpecification end

@enum Spec TemperatureSpecification BetaSpecification PressureSpecification MolarVolumeSpecification MolarEntropySpecification MolarEnthalpySpecification MolarInternalEnergySpecification

struct Specification{T}
    type::Spec
    val::T
    idx::Int
end

function (spec::Spec)(val::T)
    if spec == BetaSpecification
        throw(ArgumentError("BetaSpecification requires passing a pair of an phase index and a corresponding phase fraction (1 => 0.0)"))
    end
    return Specification(spec,val,0)
end

function (spec::Spec)(v::Pair{Int,T}) where T
    if spec != BetaSpecification
        throw(ArgumentError("invalid syntax for $(spec) . use $(spec)($(value(v))) instead"))
    end
    idx, val = v
    return Specification(spec,val,idx)
end

function r_Jacobian(eq::EqObject,x,spec::Specification)
    cache = similar(x)
    cache .= 0
    return r_Jacobian!(cache,eq,x,spec)
end

function r_Jacobian!(cache,eq,x,spec::Specification)
    cache .= 0
    sp = spec.spec
    if sp == TemperatureSpecification
        T = spec.val
        r = first(x) - T
        cache[begin] = 1
    elseif sp == BetaSpecification
        idx = spec.idx
        β = spec.val
        βi = viewlast(x,eq.Nphases)
        r = βi[idx] - β
        Ji = viewlast(cache,eq.Nphases)
        Ji[idx] = 1
    elseif sp == PressureSpecification
        p = spec.val
        r = eq.p0 - p
        cache[begin] = eq.dpdT0
        cc = @view(cache[2:eq.Ncomponents + 1])
        cc .= eq.dpdρᵢ0
    elseif sp == MolarVolumeSpecification
        v = spec.val
        T = x[begin]
        x_rest = @view(x[2:end])
        ρii = similar(x,length(eq.Nphases))
        
        for i in 1:Nphases
            ρij = viewn(x_rest,eq.Ncomponents,l)
            ρii[i] = sum(ρij)
        end
        βi = viewlast(x,eq.Nphases)
        ∑vβ = zero(eltype(x_rest))
        for i in 1:Nphases
            ∑vβ += βi[i]/ρii[i] 
        end
        Ji = viewlast(cache,eq.Nphases)
        Ji .= 1 ./ ρii
        #TODO: Jrow.segment(1+iphase*sidecar.Ncomponents, sidecar.Ncomponents) = -betas[iphase]/rho_phase[iphase]/rho_phase[iphase];
        r = v - ∑vβ
    elseif spec == MolarEntropySpecification
        T = x[begin]
    end
    return r,cache
end

struct UnpackedVariables{X}
    T::X
    ρᵢ::Vector{Vector{X}}
    β::Vector{X}
end

function pack(init::UnpackedVariables{T})
    Nphases = length(init.β)
    Ncomponents = length(init.ρᵢ[1])
    x = zeros(T,1+Nphases*(Ncomponents + 1))
    x[begin] = init.T
    x_rest = @view(x[2:end])
    for i in 1:Nphases
        xi = viewn(x_rest,Ncomponents,l)
        ρi = init.ρᵢ[i]
        xi .= ρi
    end
    x_last = viewlast(x,Nphases)
    x_last .= init.β
    return x
end

function unpack(x)

end

function residual_and_jacobian(EqObject)

end

end #module