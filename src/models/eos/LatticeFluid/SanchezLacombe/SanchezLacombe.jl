#=
struct SanchezLacombeParam <: EoSParam
    Mw::SingleParam{Float64}
    r::SingleParam{Float64}
    v_ast::SingleParam{Float64}
    omega::SingleParam{Float64}
    epsilon_ast::PairParam{Float64}
end

abstract type SanchezLacombeModel <: SAFTModel end
@newmodel SanchezLacombe SanchezLacombeModel SanchezLacombeParam
const SL = SanchezLacombe

const k_B = 1.380649e-23 # m2 kg s-2 K-1
#im supposing that this is the expression for total helmholtz energy, including ideal terms
function eos(model::SanchezLacombeModel,V,T,z=SA[1.0])
    rᵢ = model.params.r.values
    ϵᵢⱼ∗ = model.params.epsilon_ast.values
    ϵᵢ∗ = model.params.epsilon_ast.diagvalues
    ω = model.params.omega.values
    v∗ = v_mix(model,V,T,z)
    Nᵢ = z
    N = ∑(Nᵢ)
    r = dot(Nᵢ,rᵢ)/N
    ϵ∗ = eps_mix(model,V,T,z)
    T∗ = ϵ∗/k_B
    T̃ = T/T∗
    V∗ = r*N*v∗
    ṽ = V/V∗
    ṽinv = 1/ṽ
    res = (ṽ-1)*log(one(ṽ)-ṽ) + log(ṽinv)/r
    #this step could be optimized
    #also, this fails in the case of Ni[i] = 0
    ϕᵢ = rᵢ.*Nᵢ ./ (r*N)
    res += sum(log.(ϕᵢ ./ω) .* ϕᵢ ./ rᵢ) 
    res *= T̃
    res -= ṽinv
    return res
end

    
function eps_mix(model::SanchezLacombeModel,V,T,z=SA[1.0])
    rᵢ = model.params.r.values
    ϵᵢⱼ∗ = model.params.epsilon_ast.values
    ϵᵢ∗ = model.params.epsilon_ast.diagvalues
    Nᵢ = z
    N = ∑(Nᵢ)
    r = dot(Nᵢ,rᵢ)/N
    Σei = ∑(rᵢ[i]*Nᵢ[i]*ϵᵢ∗[i] for i in 1:length(z))/(r*N)
    return Σei
end

function v_mix(model::SanchezLacombeModel,V,T,z=SA[1.0])
    return 1
end

=#