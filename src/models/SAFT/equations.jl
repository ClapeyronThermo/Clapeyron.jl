function lb_volume(model::SAFTModel, z = SA[1.0])
    seg = model.params.segment.values
    σᵢᵢ = model.params.sigma.values
    val = π/6*N_A*sum(z[i]*seg[i]*σᵢᵢ[i,i]^3 for i in 1:length(z))
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
    ϵij = model.params.epsilon.values
    return view(ϵij,diagind(ϵij))
end

function p_scale(model::SAFTModel,z=SA[1.0])
    ϵ = model.params.epsilon.values
    σᵢᵢ = model.params.sigma.values
    val =  sum(z[i]*σᵢᵢ[i,i]^3/ϵ[i,i] for i in 1:length(z))*N_A/R̄
    return 1/val
end


