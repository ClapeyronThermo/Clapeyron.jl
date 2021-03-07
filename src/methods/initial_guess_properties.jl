#=
function x0_volume(model::EoS,z; phase = "unknown")
    if phase == "unknown" || is_liquid(phase)
        if typeof(model)<:SAFTgammaMie
            x0 = [log10(π/6*N_A*sum(z[i]*sum(model.group_multiplicities[i][k]*model.params.segment[k]*model.params.shapefactor[k]*model.params.sigma[k]^3 for k in @groups(i)) for i in @comps)/0.8)]
        elseif typeof(model)<:SAFT
            x0 = [log10(π/6*N_A*sum(z[i]*model.params.segment[i]*model.params.sigma[i]^3 for i in model.components)/0.8)]
        elseif typeof(model)<:Cubic
            x0 = [log10(sum(z[i]*z[j]*model.params.b[union(i,j)] for i in model.components for j in model.components)/0.8)]
        end
    elseif is_vapour(phase)
        if typeof(model)<:SAFTgammaMie
            x0 = [log10(π/6*N_A*sum(z[i]*sum(model.group_multiplicities[i][k]*model.params.segment[k]*model.params.shapefactor[k]*model.params.sigma[k]^3 for k in @groups(i)) for i in @comps)/1e-2)]
        elseif typeof(model)<:SAFT
            x0 = [log10(π/6*N_A*sum(z[i]*model.params.segment[i]*model.params.sigma[i]^3 for i in model.components)/1e-2)]
        elseif typeof(model)<:Cubic
            x0 = [log10(sum(z[i]*z[j]*model.params.b[union(i,j)] for i in model.components for j in model.components)/1e-2)]
        end
    elseif is_supercritical(phase)
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
=#

#=x0_volume=#
# function x0_volume(model::SAFTgammaMie,z; phase = "unknown")
#     if phase == "unknown" || is_liquid(phase)
#             x0 = [log10(π/6*N_A*sum(z[i]*sum(model.group_multiplicities[i][k]*model.params.segment[k]*model.params.shapefactor[k]*model.params.sigma[k]^3 for k in @groups(i)) for i in @comps)/0.8)]
#     elseif is_vapour(phase)
#             x0 = [log10(π/6*N_A*sum(z[i]*sum(model.group_multiplicities[i][k]*model.params.segment[k]*model.params.shapefactor[k]*model.params.sigma[k]^3 for k in @groups(i)) for i in @comps)/1e-2)]
#     elseif is_supercritical(phase)
#             x0 = [log10(π/6*N_A*sum(z[i]*sum(model.group_multiplicities[i][k]*model.params.segment[k]*model.params.shapefactor[k]*model.params.sigma[k]^3 for k in @groups(i)) for i in @comps)/0.5)]
#     end
#     return x0
# end

function x0_volume(model::SAFTModel,z=[1.0]; phase = "unknown")
    seg = model.params.segment.values
    σᵢᵢ = model.params.sigma.diagvalues
    val = π/6*N_A*sum(z[i]*seg[i]*σᵢᵢ[i]^3 for i in @comps)
    
    if phase == "unknown" || is_liquid(phase)
        x0val = val/0.8
    elseif is_vapour(phase)
        x0val = val/1e-2
    elseif is_supercritical(phase)
        x0val = val/0.5
    end
    return [log10(x0val)]
end

# function x0_volume(model::LJSAFT,z; phase = "unknown")
#     if phase == "unknown" || is_liquid(phase)
#         x0 = [log10(π/6*sum(z[i]*model.params.segment[i]*model.params.b[i] for i in model.components)/0.8)]
#     elseif is_vapour(phase)
#         x0 = [log10(π/6*sum(z[i]*model.params.segment[i]*model.params.b[i] for i in model.components)/1e-2)]
#     elseif is_supercritical(phase)
#         x0 = [log10(π/6*sum(z[i]*model.params.segment[i]*model.params.b[i] for i in model.components)/0.5)]
#     end
#     return x0
# end

function x0_volume(model::CubicModel,z; phase = "unknown")
    b = model.params.b.values
    val = mapreduce(+,*,b,z)
    
    if phase == "unknown" || is_liquid(phase)
        x0val = val/0.8
    elseif is_vapour(phase)
        x0val = val/1e-2
    elseif is_supercritical(phase)
        x0val = val/0.5
    end
    return [log10(x0val)]
end


#=x0_sat_pure=#

# function x0_sat_pure(model::SAFTgammaMie)
#     x0    = [log10(π/6*N_A*sum(model.group_multiplicities[model.components[1]][k]*model.params.segment[k]*model.params.shapefactor[k]*model.params.sigma[k]^3 for k in @groups(model.components[1]))/0.5),
#     log10(π/6*N_A*sum(model.group_multiplicities[model.components[1]][k]*model.params.segment[k]*model.params.shapefactor[k]*model.params.sigma[k]^3 for k in @groups(model.components[1]))/1e-3)]
# end

# function x0_sat_pure(model::SAFTVRQMie)
#     x0    = [log10(π/6*N_A*model.params.segment[model.components[1]]*model.params.sigma[model.components[1]]^3/0.2),
#     log10(π/6*N_A*model.params.segment[model.components[1]]*model.params.sigma[model.components[1]]^3/1e-3)]
# end

# function x0_sat_pure(model::LJSAFT)
#     x0    = [log10(π/6*model.params.segment[model.components[1]]*model.params.b[model.components[1]]/0.5),
#     log10(π/6*model.params.segment[model.components[1]]*model.params.b[model.components[1]]/1e-3)]
# end

function x0_sat_pure(model::SAFTModel)
    seg = only(model.params.segment.values)
    σ = only(model.params.sigma.values)
    val = π/6*N_A*seg*σ^3
    
    x0  = [val/0.5,val/1e-3]
    return log10.(x0)
end

function x0_sat_pure(model::CubicModel)
    b = only(model.params.b.values)
    x0 = [b/0.9,b/1e-4]
    return log10.(x0)
end

##=lb_volume=#
#
#lb_volume(model::SAFTgammaMie,z; phase = "unknown") = [log10(π/6*N_A*sum(z[i]*sum(model.group_multiplicities[i][k]*model.params.segment[k]*model.params.shapefactor[k]*model.params.sigma[k]^3 for k in @groups(i)) for i in @comps)/1)]
#lb_volume(model::LJSAFT,z; phase = "unknown") = [log10(π/6*sum(z[i]*model.params.segment[i]*model.params.b[i] for i in model.components)/1)]

function lb_volume(model::SAFTModel, z; phase = "unknown")
    seg = model.params.segment.values
    σᵢᵢ = model.params.sigma.diagvalues
    val = π/6*N_A*sum(z[i]*seg[i]*σᵢᵢ[i]^3 for i in @comps)
    return [log10(val)]
end
function lb_volume(model::CubicModel,z; phase = "unknown") 
    b = model.params.b.values
    val = mapreduce(+,*,b,z)
    return [log10(val)]

end

lb_volume(model::IAPWS95, z; phase = "unknown") = [-5.0]

#=scale_sat_pure=#

# function scale_sat_pure(model::SAFTgammaMie)
#     m̄  = sum(model.group_multiplicities[model.components[1]][k]*model.params.segment[k]*model.params.shapefactor[k] for k in @groups(model.components[1]))
#         σ̄3 = sum(model.group_multiplicities[model.components[1]][k]*model.group_multiplicities[model.components[1]][l]*
#                  model.params.segment[k]*model.params.segment[l]*
#                  model.params.shapefactor[k]*model.params.shapefactor[l]*
#                  model.params.sigma[union(k,l)]^3 for k in @groups(model.components[1]) for l in @groups(model.components[1]))/m̄^2
#         ϵ̄  = T_crit_pure(model)
#         p_scale    = σ̄3*N_A/R̄/ϵ̄
#         μ_scale    = 1/R̄/ϵ̄
#     return p_scale,μ_scale
# end

# function scale_sat_pure(model::LJSAFT)
#     p_scale    = model.params.b[model.components[1]]/N_A/R̄/model.params.T[model.components[1]]
#     μ_scale    = 1/R̄/model.params.T[model.components[1]]
#     return p_scale,μ_scale
# end

function scale_sat_pure(model::SAFTModel)
    ϵ = only(model.params.epsilon.values[1])
    σ = only(model.params.sigma.values)

    p_scale    = σ^3*N_A/R̄/ϵ
    μ_scale    = 1/R̄/ϵ
    return p_scale,μ_scale
end

function scale_sat_pure(model::CubicModel)
    a = only(model.params.a.values)
    b = only(model.params.b.values)
    p_scale    = b^2/a*27
    μ_scale    = 27*b/8/a
    return p_scale,μ_scale
end

#=x0_crit_pure=#

# x0_crit_pure(model::SAFTgammaMie) = [2, log10(π/6*N_A*sum(model.group_multiplicities[model.components[1]][k]*model.params.segment[k]*model.params.shapefactor[k]*model.params.sigma[k]^3 for k in @groups(model.components[1]))/0.3)]
# x0_crit_pure(model::LJSAFT) = [1.5, log10(π/6*model.params.segment[model.components[1]]*model.params.b[model.components[1]]/0.3)]

x0_crit_pure(model::SAFTModel) = [1.5, log10(π/6*N_A*model.params.segment.values[1]*model.params.sigma.values[1]^3/0.3)]
# x0_crit_pure(model::Cubic) = [1.0, log10(model.params.b[model.components[1]]/0.3)]
# x0_crit_pure(model::CPA) = [2, log10(model.params.b[model.components[1]]/0.3)]



#=T_crit_pure=#
# function T_crit_pure(model::SAFTgammaMie)
#     m̄ = sum(model.group_multiplicities[model.components[1]][k]*model.params.segment[k]*model.params.shapefactor[k] for k in @groups(model.components[1]))
#     return sum(model.group_multiplicities[model.components[1]][k]*model.group_multiplicities[model.components[1]][l]*
#                model.params.segment[k]*model.params.segment[l]*
#                model.params.shapefactor[k]*model.params.shapefactor[l]*
#                model.params.epsilon[union(k,l)] for k in @groups(model.components[1]) for l in @groups(model.components[1]))/m̄^2
# end
# T_crit_pure(model::LJSAFT) = model.params.T[model.components[1]]

T_crit_pure(model::SAFTModel) = model.params.epsilon.values[1]

# T_crit_pure(model::Cubic) = model.params.a[model.components[1]]/model.params.b[model.components[1]]/8.314*8/27
