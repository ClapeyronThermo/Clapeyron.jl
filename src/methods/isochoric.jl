#isochoric formalism
#=
ρ = 1/V
xi = ρi/ρ
Ψ(model,T,ρi) = ρ*a_res(model,V,T,z)

#we also have the following relation:
ρ = sum(z)/V

we can transform from Clapeyron formalism to isochoric formalism by using V = 1:
we suppose ∑z = ρ, then:
zi = xi .* ∑z = ρi/ρ * ρ = ρi
V = sum(z)/ρ = ρ/ρ = 1
=#

function Psi_ideal(model,T,ρᵢ)
    ρ = sum(ρᵢ)
    V = one(ρ)
    return Rgas(model)*T*ρ*a_ideal(model,one(ρ),T,ρᵢ)  + reference_state_eval(model,V,T,ρᵢ)
end

const Ψ_ideal = Psi_ideal

function Psi(model,T,ρᵢ)
    ρ = sum(ρᵢ)
    V = one(ρ)
    return Rgas(model)*T*ρ*a_eos(model,V,T,ρᵢ) + reference_state_eval(model,V,T,ρᵢ)
end

function Psi(model::IdealModel,T,ρᵢ)
    ρ = sum(ρᵢ)
    V = one(ρ)
    return Rgas(model)*T*ρ*a_ideal(model,V,T,ρᵢ) + reference_state_eval(model,V,T,ρᵢ)
end

const Ψ_eos = Psi

function Psi_res(model,T,ρᵢ)
    ρ = sum(ρᵢ)
    return Rgas(model)*T*ρ*a_res(model,one(ρ),T,ρᵢ)
end

const Ψ_res = Psi_res

#returns gradient(Ψ_eos) and hessian(Ψ_eos)
function Ψ_grad_and_hessian_res(model,T,ρᵢ)
    result = DiffResults.HessianResult(ρᵢ)
    return Ψ_grad_and_hessian_res!(result,model,T,ρᵢ)
end

function Ψ_grad_and_hessian_res!(result,model,T,ρᵢ)
    Ψ(ρ) = Psi_res(model,T,ρ)
    result = ForwardDiff.hessian!(result,Ψ,ρᵢ)
    ∇Ψ = DiffResults.gradient(result)
    HΨ = DiffResults.hessian(result)
    return ∇Ψ, HΨ
end

function Ψ_grad_and_hessian(model,T,ρᵢ)
    result = DiffResults.HessianResult(ρᵢ)
    return Ψ_grad_and_hessian!(result,model,T,ρᵢ)
end

function Ψ_grad_and_hessian!(result,model,T,ρᵢ)
    Ψ(ρ) = Psi(model,T,ρ)
    ∇Ψ, HΨ = Ψ_grad_and_hessian_res!(result,model,T,ρᵢ)
    return ∇Ψ, HΨ
end

function Ψ_hessian_res(model,T,ρᵢ)
    Ψ(ρ) = Psi_res(model,T,ρ)
    HΨ = ForwardDiff.hessian(Ψ,ρᵢ)
    return HΨ
end

function Ψ_hessian(model,T,ρᵢ)
    HΨr = Ψ_hessian_res(model,T,ρᵢ)
    RT = Rgas(model)*T
    for i in 1:length(ρᵢ)
        HΨr[i,i] += RT/ρᵢ[i]
    end
    return HΨr
end

function Ψ_grad(model,T,ρᵢ)
    TT = gradient_type(model,T,ρᵢ)
    similar(TT,length(ρᵢ))
    return Ψ_grad!(similar(TT,length(ρᵢ)),model,T,ρᵢ)
end

function Ψ_grad!(∇Ψ,model,T,ρᵢ)
    Ψ(ρ) = Psi(model,T,ρ)
    ForwardDiff.gradient!(∇Ψ,Ψ,dΨrdT)
    return ∇Ψ
end

function dΨdT_grad_res(model,T,ρᵢ)
    TT = gradient_type(model,T,ρᵢ)
    similar(TT,length(ρᵢ))
    return dΨdT_grad_res!(similar(TT,length(ρᵢ)),model,T,ρᵢ)
end

function dΨdT_grad_res!(∇Ψ,model,T,ρᵢ)
    dΨrdT(ρ) = last(Solvers.f∂f(_T -> Psi_res(model,_T,ρ),T))
    ForwardDiff.gradient!(∇Ψ,Ψ,dΨrdT)
    return ∇Ψ
end

function d2ΨdT2_grad_res(model,T,ρᵢ)
    TT = gradient_type(model,T,ρᵢ)
    similar(TT,length(ρᵢ))
    return d2ΨdT2_grad_res!(similar(TT,length(ρᵢ)),model,T,ρᵢ)
end

function d2ΨdT2_grad_res!(∇Ψ,model,T,ρᵢ)
    d2ΨrdT2(ρ) = last(Solvers.f∂f∂2f(_T -> Psi_res(model,_T,ρ),T))
    ForwardDiff.gradient!(∇Ψ,Ψ,d2ΨrdT2)
    return ∇Ψ
end 

function dΨdT_grad(model,T,ρᵢ)
    TT = gradient_type(model,T,ρᵢ)
    similar(TT,length(ρᵢ))
    return dΨdT_grad!(similar(TT,length(ρᵢ)),model,T,ρᵢ)
end

function dΨdT_grad!(∇Ψ,model,T,ρᵢ)
    dΨrdT(ρ) = last(Solvers.f∂f(_T -> Psi(model,_T,ρ),T))
    ForwardDiff.gradient!(∇Ψ,Ψ,dΨrdT)
    return ∇Ψ
end

function d2ΨdT2_grad(model,T,ρᵢ)
    TT = gradient_type(model,T,ρᵢ)
    similar(TT,length(ρᵢ))
    return d2ΨdT2_grad!(similar(TT,length(ρᵢ)),model,T,ρᵢ)
end

function d2ΨdT2_grad!(∇Ψ,model,T,ρᵢ)
    d2ΨrdT2(ρ) = last(Solvers.f∂f∂2f(_T -> Psi(model,_T,ρ),T))
    ForwardDiff.gradient!(∇Ψ,Ψ,d2ΨrdT2)
    return ∇Ψ
end 
