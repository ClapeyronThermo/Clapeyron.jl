using Serialization

const FeANNMieConsts::@NamedTuple{v_kernel_0::Vector{Float64}, bias_0::Vector{Float64}, v_kernel_1::Vector{Float64}, bias_1::Vector{Float64}, v_kernel_2::Vector{Float64}, bias_2::Vector{Float64}, v_kernel_3::Vector{Float64}, bias_3::Vector{Float64}, kernel_helmholtz::Vector{Float64}} = Serialization.deserialize("mie_feann.jls")

struct FeANNMie{T} <: EoSModel
    lambda_r::T
    lambda_a::T
end

function a_res(model::FeANNMie,V,T,z)
    ρ⃰ = sum(z)/V
    T⃰ = T
    λa, λr = model.lambda_a,model.lambda_r
    C = Cλ_mie(λa, λr)
    α = C*(1/(λa - 3) - 1/(λr - 3))
    x1 = zeros(Base.promote_eltype(model,V,T,z),45)
    x2 = similar(x2)
    cache = (x1,x2)
    a = a_fe_ann(α, T⃰, ρ⃰,cache)
    a0 = a_fe_ann(α, T⃰, zero(ρ⃰), cache)
    return a - a0
end

function a_fe_ann(α, T⃰, ρ⃰,cache)
    # Input row vector: 1 x 3
    x1,x2,x3 = promote(α, ρ⃰, 1.0 / T⃰)
    x = SVector(x1,x2,x3)

    mul_cache,z = cache
    # Layer 0 
    kernel_0 = resize(FeANNMieConsts.v_kernel_0,(3,45))
    bias_0 = FeANNMieConsts.bias_0
    mul!(mul_cache,kernel_0,x)
    z .= tanh.(mul_cache .+ bias_0)
    # Layer 1
    kernel_1 = resize(FeANNMieConsts.v_kernel_1,(45,45))
    bias_1 = FeANNMieConsts.bias_1
    mul!(mul_cache,kernel_1,z)
    z .= tanh.(mul_cache .+ bias_1)
    # Layer 2
    kernel_2 = resize(FeANNMieConsts.v_kernel_2,(45,45))
    bias_2 = FeANNMieConsts.bias_2
    mul!(mul_cache,kernel_2,z)
    z .= tanh.(mul_cache .+ bias_2)
    # Layer 3
    kernel_3 = resize(FeANNMieConsts.v_kernel_3,(45,45))
    bias_3 = FeANNMieConsts.bias_3
    mul!(mul_cache,kernel_3,z)
    z .= tanh.(mul_cache .+ bias_3)

    return dot(z,FeANNMieConsts.kernel_helmholtz)/T⃰
end
