function lb_volume(model::SAFTModel, z = SA[1.0])
    seg = model.params.segment.values
    σ = model.params.sigma.values
    val = π/6*N_A*sum(z[i]*seg[i]*σ[i,i]^3 for i in 1:length(z))
    return val
end

function x0_crit_pure(model::SAFTModel)
    lb_v = lb_volume(model)
    (2.0, log10(lb_v/0.3))
end

function T_scale(model::SAFTModel,z=SA[1.0])
    ϵ = model.params.epsilon.values
    return prod(ϵ[i,i]^z[i] for i in 1:length(z))^(1/sum(z))
end

function T_scales(model::SAFTModel)
    ϵ =diagvalues(model.params.epsilon)
end

function p_scale(model::SAFTModel,z=SA[1.0])
    ϵ = model.params.epsilon.values
    σ = model.params.sigma.values
    val =  sum(z[i]*σ[i,i]^3/ϵ[i,i] for i in 1:length(z))*N_A/R̄
    return 1/val
end

function packing_fraction(model::SAFTModel,V,T,z)
    σ = diagvalues(model.params.sigma)
    m = model.params.segment.values
    x = z ./ ∑(z)
    return π*N_A/V/6*∑(x.*m.*σ.^3)
end

function antoine_coef(model::SAFTModel)
    m = model.params.segment.values[1]
    A = 2.3461144513376593+0.27679968565666935*m
    B = exp(1.7330494260220226 + 0.6185684341246401*log(m))
    C = 0.018524160155803788 - 0.19222021003570597*log(m)
    return A,B,C
end    

## Association overloads required to support association

@inline function assoc_similar(model::Union{SAFTModel,CPAModel},::Type{𝕋}) where 𝕋
    assoc_similar(model.params.bondvol.values,𝕋)
end

#recombine! utilities
function recombine_saft!(model::SAFTModel)
    sigma = model.params.sigma
    epsilon = model.params.epsilon
    sigma = sigma_LorentzBerthelot!(sigma)
    epsilon = epsilon_LorentzBerthelot!(epsilon)
    return model
end