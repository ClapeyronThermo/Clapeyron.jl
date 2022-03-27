function lb_volume(model::SAFTModel, z = SA[1.0])
    seg = model.params.segment.values
    σᵢᵢ = model.params.sigma.diagvalues
    val = π/6*N_A*sum(z[i]*seg[i]*σᵢᵢ[i]^3 for i in 1:length(z))
    return val
end

function x0_crit_pure(model::SAFTModel)
    lb_v = lb_volume(model)
    (2.0, log10(lb_v/0.3))
end

function T_scale(model::SAFTModel,z=SA[1.0])
    ϵ = model.params.epsilon.diagvalues
    return prod(ϵ[i]^z[i] for i in 1:length(z))^(1/sum(z))
end

function T_scales(model::SAFTModel)
    ϵ = model.params.epsilon.diagvalues
end

function p_scale(model::SAFTModel,z=SA[1.0])
    ϵ = model.params.epsilon.diagvalues
    σᵢᵢ = model.params.sigma.diagvalues
    val =  sum(z[i]*σᵢᵢ[i]^3/ϵ[i] for i in 1:length(z))*N_A/R̄
    return 1/val
end


