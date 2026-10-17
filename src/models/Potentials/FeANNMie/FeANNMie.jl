#=most of this implementation is taken from https://github.com/gustavochm/Square-Peg-Mie

Weights were read by MsgPack:

data_dict =open("Clapeyron/src/models/Potentials/FeANNMie/checkpoint","r") do io
          data = MsgPack.unpack(io)
       end

params = data_dict["params]

function unpack_flax_tensor(ext::MsgPack.Extension)
    @assert ext.type == 1 "unexpected extension type code"
    io = IOBuffer(ext.data)
    shape, dtype, raw = MsgPack.unpack(io)
    T = dtype == "float32" ? Float32 :
        dtype == "float64" ? Float64 :
        error("unsupported dtype $dtype")
    A = reinterpret(T, raw)
    return reshape(collect(A), Tuple(shape)...)
end

function save_compatible(path::AbstractString, data; target_version=UInt8(24))
    buf = IOBuffer()
    serialize(buf, data)
    bytes = take!(buf)
    bytes[4] = target_version   # 1-based index 4 == header offset 3
    write(path, bytes)
end

b0 = unpack_flax_tensor(params["hidden_layers_0"]["bias"])
b1 = unpack_flax_tensor(params["hidden_layers_1"]["bias"])
b2 = unpack_flax_tensor(params["hidden_layers_2"]["bias"])
W0 = unpack_flax_tensor(params["hidden_layers_0"]["kernel"])
W1 = unpack_flax_tensor(params["hidden_layers_1"]["kernel"])
W2 = unpack_flax_tensor(params["hidden_layers_2"]["kernel"])
Wh = unpack_flax_tensor(params["helmholtz_layer"]["kernel"])

nt_data = (;W0,b0,W1,b1,W2,b2,Wh)
save_compatible("Clapeyron/src/models/Potentials/FeANNMie/params_anns.jls",nt_data)
=#
using Serialization

function load_fe_ann_mie_consts(solid::Bool)
    if solid
        p = open(deserialize, joinpath(@__DIR__, "params_anns.jls"))
    else
        p = open(deserialize, joinpath(@__DIR__, "params_ann.jls"))
    end
    #files extracted from https://github.com/gustavochm/Square-Peg-Mie/tree/main/ann_models/feann_eos
    #we do additional reshaping here to make it work with julia conventions.
    return (
        solid = solid,
        W0t = reshape(p.W0, 30,  3),
        b0  = p.b0,
        W1t = reshape(p.W1, 30, 30),
        b1  = p.b1,
        W2t = reshape(p.W2, 30, 30),
        b2  = p.b2,
        Wh  = vec(p.Wh),
    )
end

struct FeANNMie{T,S} <: EoSModel
    lambda_r::T
    lambda_a::T
    weights::S
end

is_splittable(::FeANNMie) = false

"""

    FeANNMie(λr, λa; solid = false)


## Input arguments

- `λa`: Real - Atractive range parameter (no units)
- `λr`: Real - Repulsive range parameter (no units)

## Keyword arguments
- `solid`: if set to `true`, it will load the weights for the `Fe-ANN(s)-Mie` model instead of the weights for the `Fe-ANN-Mie`.

## Description

Free Energy Artificial Neural Network equation of state for the Mie fluid.
There are two available models, that can be selected with the `solid` keyword argument:
- Fe-ANN-Mie: model trained on VLE data of the Mie fluid
- Fe-ANN(S)-Mie: Trained model using fluid and crystal Mie data

## References

1. Chaparro, G., & Müller, E. A. (2023). Development of thermodynamically consistent machine-learning equations of state: Application to the Mie fluid. The Journal of Chemical Physics, 158(18), 184505. [doi:10.1063/5.0146634](https://doi.org/10.1063/5.0146634)
1. Chaparro, G., & Müller, E. A. (2024). Development of a Helmholtz free energy equation of state for fluid and solid phases via artificial neural networks. Communications Physics, 7(1). [doi:10.1038/s42005-024-01892-3]https://doi.org/10.1038/s42005-024-01892-3)
1. Chaparro, G., & Müller, E. A. (2025). Fitting a square peg in a round hole: parameterisation of quasi-spherical molecules employing the Mie potential. Molecular Systems Design & Engineering, 10(8), 620–634. [doi:10.1039/d5me00048c](https://doi.org/10.1039/d5me00048c)
"""
function FeANNMie(λr, λa; solid::Bool = false)
    λr, λa = promote(λr, λa)
    weights = load_fe_ann_mie_consts(solid)
    return FeANNMie(λr, λa, weights)
end

function Base.show(io::IO,model::FeANNMie)
    print(io,"FeANNMie(λr = $(model.lambda_r), λa = $(model.lambda_a)")
    model.weights.feanns && print(io,", solid = true")
    print(io,")")
end

function Base.show(io::IO,::MIME"text/plain",model::FeANNMie)
    print(io,"FeANNMie(λr = $(model.lambda_r), λa = $(model.lambda_a)")
    model.weights.solid && print(io,", solid = true")
    print(io,")")
end

function mie_α_const(λa, λr)
    C = Cλ_mie(λa, λr)
    α = C*(1/(λa - 3) - 1/(λr - 3))
    return α
end

mie_α_const(model::FeANNMie) = mie_α_const(model.lambda_a,model.lambda_r)
function a_res(model::FeANNMie,V,T,z)
    ρ⃰ = sum(z)/V
    T⃰ = T
    λa, λr = model.lambda_a,model.lambda_r
    α = mie_α_const(λa, λr)
    x1 = zeros(Base.promote_eltype(model,V,T,z),30)
    x2 = similar(x1)
    cache = (x1,x2)
    weights = model.weights
    a = a_fe_ann(weights, α, T⃰, ρ⃰,cache)
    a0 = a_fe_ann(weights, α, T⃰, zero(ρ⃰), cache)
    return a - a0
end

function a_fe_ann(weights, α, T⃰, ρ⃰)
    x1 = zeros(Base.promote_eltype(α, T⃰, ρ⃰),30)
    x2 = similar(x1)
    return a_fe_ann(weights, α, T⃰, ρ⃰,(x1,x2))
end

function a_fe_ann(weights, α, T⃰, ρ⃰, cache)
    mul_cache, z = cache
    x = SVector(promote(α, ρ⃰, 1.0 / T⃰))

    mul!(mul_cache, weights.W0t, x)
    z .= tanh.(mul_cache .+ weights.b0)

    mul!(mul_cache, weights.W1t, z)
    z .= tanh.(mul_cache .+ weights.b1)

    mul!(mul_cache, weights.W2t, z)
    z .= tanh.(mul_cache .+ weights.b2)

    return dot(z, weights.Wh) / T⃰
end

"""
    galliero_tc(ms, α)

Calculate the critical temperature using the Galliero's CS correlation.
Requires `ms` (chain length) and `α` (van der Waals constant for the mie fluid) 

## References
1. Hoang, H., Delage-Santacreu, S., & Galliero, G. (2017). Simultaneous description of equilibrium, interfacial, and transport properties of fluids using a Mie chain coarse-grained force field. Industrial & Engineering Chemistry Research, 56(32), 9213–9226. doi:10.1021/acs.iecr.7b01397
"""
function galliero_tc(ms, α)

    # Table 1
    A0, A1, A2, A3, A4, A5, A6, A7, A8 = -4.0219, -0.6969, 3.1689, -0.5822,  -0.3909, 1.4849, 0.7195, -3.1836, 0.2520

    # Eq. 10 
    sqrt_ms = sqrt(ms)
    Tc = (A0 * exp(A1 * sqrt_ms) + A2) * α^(A3 * exp(A4 * ms) + A5)
    Tc += (A6 * exp(A7 / sqrt_ms) + A8)
    return Tc
end

"""
    galliero_rholsat(ms, α)

Calculate the an estimate of the saturated liquid density at Tr = 0.7 using the Galliero's CS correlation.
Requires `ms` (chain length) and `α` (van der Waals constant for the mie fluid)

## References
1. Hoang, H., Delage-Santacreu, S., & Galliero, G. (2017). Simultaneous description of equilibrium, interfacial, and transport properties of fluids using a Mie chain coarse-grained force field. Industrial & Engineering Chemistry Research, 56(32), 9213–9226. doi:10.1021/acs.iecr.7b01397
"""
function galliero_rholsat(ms, α)

    # Table 1
    B0, B1, B2, B3, B4, B5, B6, B7, B8 = 0.3198, -0.6815, 0.2532, -1.9621, -0.4389, -2.5033, 0.1568, -0.1756, 0.5972

    # Eq. 11
    rholsat = (B0 * exp(B1 * ms) + B2) * exp((B3 * exp(B4 * ms) + B5) * α^1.75) 
    rholsat += (B6 * exp(B7 * ms) + B8)
    rholsat /= ms
    return rholsat
end

"""
    galliero_acentric_factor(ms, α)

Calculate the an estimate of acentric factor using the Galliero's CS correlation.
Requires `ms` (chain length) and `α` (van der Waals constant for the mie fluid) 

## References
1. Hoang, H., Delage-Santacreu, S., & Galliero, G. (2017). Simultaneous description of equilibrium, interfacial, and transport properties of fluids using a Mie chain coarse-grained force field. Industrial & Engineering Chemistry Research, 56(32), 9213–9226. doi:10.1021/acs.iecr.7b01397
"""
function galliero_acentric_factor(ms, α)

    # Table 1
    C0, C1, C2, C3 = -0.1574, -0.0008, 0.3158, -0.4292

    # Eq. 12
    α∞ = 1/3
    sqrt_ms = sqrt(ms)
    ω = (C0 * sqrt_ms + C1) * log(α - α∞)
    ω += (C2 * sqrt_ms + C3)
    return ω
end

T_scale(::FeANNMie,z) = 1.0
Base.length(::FeANNMie) = 1
paramtype(::FeANNMie{T}) where T = T
paramtype(::Type{FeANNMie{T}}) where T = T
component_list(model::FeANNMie) = ["Mie ($(model.lambda_r)-$(model.lambda_a)) fluid"]
lb_volume(::FeANNMie,T,z) = 1/1.25
Rgas(::FeANNMie) = 1.0
idealmodel(::FeANNMie) = BasicIdeal()

function x0_volume_liquid(model::FeANNMie,p,T,z)
    α = mie_α_const(model)
    tc = galliero_tc(1.0, α)
    _1 = one(Base.promote_eltype(p,T,z))
    rhol = galliero_rholsat(1.0, α)*_1
    v0 = sum(z)/rhol
    T > 0.7*tc && return v0
    return volume_chill(model,p,T,z,v0,0.7tc)
end

function x0_crit_pure(model::FeANNMie,z)
    α = mie_α_const(model)
    Tr = galliero_tc(1.0, α) #T_scale is already one
    return Tr, log10(1/galliero_rholsat(1.0, α))
end

export FeANNMie