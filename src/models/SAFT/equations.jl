function lb_volume(model::SAFTModel, z = SA[1.0])
    m = model.params.segment.values
    σ = model.params.sigma.values
    val = π/6*N_A*sum(z[i]*m[i]*σ[i,i]^3 for i in 1:length(z))
    return val
end

function x0_crit_pure(model::SAFTModel)
    lb_v = lb_volume(model)
    (2.0, log10(lb_v/0.3))
end

function saft_lorentz_berthelot(params)
    k = get(params,"k",nothing)
    l = get(params,"l",nothing)
    sigma,epsilon = params["sigma"],params["epsilon"]
    params["sigma"] = sigma_LorentzBerthelot(sigma, l)
    params["epsilon"] = epsilon_LorentzBerthelot(epsilon, k)
    return params
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
    val = sum(z[i]*σ[i,i]^3/ϵ[i,i] for i in 1:length(z))*N_A/R̄
    return 1/val
end

function antoine_coef(model::SAFTModel)
    m = model.params.segment.values[1]
    A = 2.3461144513376593+0.27679968565666935*m
    B = exp(1.7330494260220226 + 0.6185684341246401*log(m))
    C = 0.018524160155803788 - 0.19222021003570597*log(m)
    return A,B,C
end    

## Association overloads required to support association

@inline function assoc_similar(model::EoSModel,::Type{𝕋}) where 𝕋
    assoc_similar(model.params.bondvol.values,𝕋)
end

#recombine! utilities
function recombine_saft!(model::SAFTModel,k = nothing,l = nothing)
    sigma = model.params.sigma
    epsilon = model.params.epsilon
    sigma = sigma_LorentzBerthelot!(sigma,l)
    epsilon = epsilon_LorentzBerthelot!(epsilon,k)
    return model
end