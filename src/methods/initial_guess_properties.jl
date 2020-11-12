function x0_volume(model::EoS,z; phase = "unknown")
    if phase == "unknown" || phase == "liquid"
        if typeof(model)<:SAFTgammaMie
            x0 = [log10(π/6*N_A*sum(z[i]*sum(model.group_multiplicities[i][k]*model.params.segment[k]*model.params.shapefactor[k]*model.params.sigma[k]^3 for k in @groups(i)) for i in @comps)/0.8)]
        elseif typeof(model)<:SAFT
            x0 = [log10(π/6*N_A*sum(z[i]*model.params.segment[i]*model.params.sigma[i]^3 for i in model.components)/0.8)]
        elseif typeof(model)<:Cubic
            x0 = [log10(sum(z[i]*z[j]*model.params.b[union(i,j)] for i in model.components for j in model.components)/0.8)]
        end
    elseif phase == "vapour"
        if typeof(model)<:SAFTgammaMie
            x0 = [log10(π/6*N_A*sum(z[i]*sum(model.group_multiplicities[i][k]*model.params.segment[k]*model.params.shapefactor[k]*model.params.sigma[k]^3 for k in @groups(i)) for i in @comps)/1e-2)]
        elseif typeof(model)<:SAFT
            x0 = [log10(π/6*N_A*sum(z[i]*model.params.segment[i]*model.params.sigma[i]^3 for i in model.components)/1e-2)]
        elseif typeof(model)<:Cubic
            x0 = [log10(sum(z[i]*z[j]*model.params.b[union(i,j)] for i in model.components for j in model.components)/1e-2)]
        end
    elseif phase == "supercritical"
        if typeof(model)<:SAFTgammaMie
            x0 = [log10(π/6*N_A*sum(z[i]*sum(model.group_multiplicities[i][k]*model.params.segment[k]*model.params.shapefactor[k]*model.params.sigma[k]^3 for k in @groups(i)) for i in @comps)/0.5)]
        elseif typeof(model)<:SAFT
            x0 = [log10(π/6*N_A*sum(z[i]*model.params.segment[i]*model.params.sigma[i]^3 for i in model.components)/0.5)]
        elseif typeof(model)<:Cubic
            x0 = [log10(sum(z[i]*z[j]*model.params.b[union(i,j)] for i in model.components for j in model.components)/0.5)]
        end
    end

    return x0
end

function lb_volume(model::EoS,z; phase = "unknown")
        if typeof(model)<:SAFTgammaMie
            lb = [log10(π/6*N_A*sum(z[i]*sum(model.group_multiplicities[i][k]*model.params.segment[k]*model.params.shapefactor[k]*model.params.sigma[k]^3 for k in @groups(i)) for i in @comps)/1)]
        elseif typeof(model)<:SAFT
            lb = [log10(π/6*N_A*sum(z[i]*model.params.segment[i]*model.params.sigma[i]^3 for i in model.components)/1)]
        elseif typeof(model)<:Cubic
            lb = [log10(sum(z[i]*z[j]*model.params.b[union(i,j)] for i in model.components for j in model.components))]
        end
    return lb
end

function x0_sat_pure(model::EoS)
    if typeof(model)<:SAFTgammaMie
        x0    = [log10(π/6*N_A*sum(model.group_multiplicities[model.components[1]][k]*model.params.segment[k]*model.params.shapefactor[k]*model.params.sigma[k]^3 for k in @groups(model.components[1]))/0.45),
                 log10(π/6*N_A*sum(model.group_multiplicities[model.components[1]][k]*model.params.segment[k]*model.params.shapefactor[k]*model.params.sigma[k]^3 for k in @groups(model.components[1]))/1e-3)]
    elseif typeof(model)<:SAFT
        x0    = [log10(π/6*N_A*model.params.segment[model.components[1]]*model.params.sigma[model.components[1]]^3/0.45),
                 log10(π/6*N_A*model.params.segment[model.components[1]]*model.params.sigma[model.components[1]]^3/1e-3)]
    elseif typeof(model)<:Cubic
        x0    = [log10(model.params.b[model.components[1]]/0.9),
                 log10(model.params.b[model.components[1]]/1e-3)]
    end
    return x0
end


function scale_sat_pure(model)
    if typeof(model)<:SAFTgammaMie
        m̄  = sum(model.group_multiplicities[model.components[1]][k]*model.params.segment[k]*model.params.shapefactor[k] for k in @groups(model.components[1]))
        σ̄3 = sum(model.group_multiplicities[model.components[1]][k]*model.group_multiplicities[model.components[1]][l]*
                 model.params.segment[k]*model.params.segment[l]*
                 model.params.shapefactor[k]*model.params.shapefactor[l]*
                 model.params.sigma[union(k,l)]^3 for k in @groups(model.components[1]) for l in @groups(model.components[1]))/m̄^2
        ϵ̄  = T_crit_pure(model)
        p_scale    = σ̄3*N_A/R̄/ϵ̄
        μ_scale    = 1/R̄/ϵ̄
    elseif typeof(model)<:SAFT
        p_scale    = model.params.sigma[model.components[1]]^3*N_A/R̄/model.params.epsilon[model.components[1]]
        μ_scale    = 1/R̄/model.params.epsilon[model.components[1]]
    elseif typeof(model)<:Cubic
        p_scale    = model.params.b[model.components[1]]^2/model.params.a[model.components[1]]*27
        μ_scale    = 27*model.params.b[model.components[1]]/8/model.params.a[model.components[1]]
    end
    return p_scale,μ_scale
end

function x0_crit_pure(model::EoS)
    if typeof(model)<:SAFTgammaMie
        x0 = [1.5, log10(π/6*N_A*sum(model.group_multiplicities[model.components[1]][k]*model.params.segment[k]*model.params.shapefactor[k]*model.params.sigma[k]^3 for k in @groups(model.components[1]))/0.3)]
    elseif typeof(model)<:SAFT
        x0 = [1.5, log10(π/6*N_A*model.params.segment[model.components[1]]*model.params.sigma[model.components[1]]^3/0.3)]
    elseif typeof(model)<:Cubic
        x0 = [1, log10(model.params.b[model.components[1]]/0.3)]
    end
    return x0
end

function T_crit_pure(model::EoS)
    if typeof(model)<:SAFTgammaMie
        m̄ = sum(model.group_multiplicities[model.components[1]][k]*model.params.segment[k]*model.params.shapefactor[k] for k in @groups(model.components[1]))
        return sum(model.group_multiplicities[model.components[1]][k]*model.group_multiplicities[model.components[1]][l]*
                   model.params.segment[k]*model.params.segment[l]*
                   model.params.shapefactor[k]*model.params.shapefactor[l]*
                   model.params.epsilon[union(k,l)] for k in @groups(model.components[1]) for l in @groups(model.components[1]))/m̄^2
    elseif typeof(model)<:SAFT
        return model.params.epsilon[model.components[1]]
    elseif typeof(model)<:Cubic
        return model.params.a[model.components[1]]/model.params.b[model.components[1]]/8.314*8/27
    end
end
