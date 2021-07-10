#=
function x0_volume(model::EoS,z; phase = :unknown)
    if phase == :unknown || is_liquid(phase)
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
# function x0_volume(model::SAFTgammaMie,z; phase = :unknown)
#     if phase == :unknown || is_liquid(phase)
#             x0 = [log10(π/6*N_A*sum(z[i]*sum(model.group_multiplicities[i][k]*model.params.segment[k]*model.params.shapefactor[k]*model.params.sigma[k]^3 for k in @groups(i)) for i in @comps)/0.8)]
#     elseif is_vapour(phase)
#             x0 = [log10(π/6*N_A*sum(z[i]*sum(model.group_multiplicities[i][k]*model.params.segment[k]*model.params.shapefactor[k]*model.params.sigma[k]^3 for k in @groups(i)) for i in @comps)/1e-2)]
#     elseif is_supercritical(phase)
#             x0 = [/0.5)]
#     end
#     return x0
# end


# function x0_volume(model::LJSAFT,z; phase = :unknown)
#     if phase == :unknown || is_liquid(phase)
#         x0 = [log10(π/6*sum(z[i]*model.params.segment[i]*model.params.b[i] for i in model.components)/0.8)]
#     elseif is_vapour(phase)
#         x0 = [log10(π/6*sum(z[i]*model.params.segment[i]*model.params.b[i] for i in model.components)/1e-2)]
#     elseif is_supercritical(phase)
#         x0 = [/0.5)]
#     end
#     return x0
# end

#=
lb_volume:

SAFTgammaMie:  log10(π/6*N_A*sum(z[i]*sum(model.group_multiplicities[i][k]*model.params.segment[k]*model.params.shapefactor[k]*model.params.sigma[k]^3 for k in @groups(i)) for i in @comps)
LJSAFT: log10(π/6*sum(z[i]*model.params.segment[i]*model.params.b[i] for i in @comps)
=#


function x0_volume_liquid(model,T,z)
    v_lb = lb_volume(model,z)
    return v_lb/0.8
end

function x0_volume_liquid(model::SAFTVRMieModel,T,z)
    v_lb = lb_volume(model,z)
    return v_lb*1.5
end

function x0_volume_gas(model,p,T,z)
    return volume_virial(model,p,T,z)
end

function x0_volume_sc(model,p,T,z)
    v_sc = lb_volume(model,z)
    return v_sc*2
end

function x0_volume(model::EoSModel,p,T,z; phase = :unknown)
    phase = Symbol(phase)
    if phase === :unknown || is_liquid(phase)
        return x0_volume_liquid(model,T,z)
    elseif is_vapour(phase)
        return x0_volume_gas(model,p,T,z)
    elseif is_supercritical(phase)
        return x0_volume_sc(model,p,T,z)
    else
        error("unreachable state on x0_volume")
    end
end


#=x0_sat_pure=#

# function x0_sat_pure(model::SAFTVRQMie)
#     x0    = [log10(π/6*N_A*model.params.segment[model.components[1]]*model.params.sigma[model.components[1]]^3/0.2),
#     log10(π/6*N_A*model.params.segment[model.components[1]]*model.params.sigma[model.components[1]]^3/1e-3)]
# end

# function x0_sat_pure(model::LJSAFT)
#     x0    = [log10(π/6*model.params.segment[model.components[1]]*model.params.b[model.components[1]]/0.5),
#     log10(π/6*model.params.segment[model.components[1]]*model.params.b[model.components[1]]/1e-3)]
# end




##=lb_volume=#
#
#lb_volume(model::LJSAFT,z; phase = :unknown) = [log10(π/6*sum(z[i]*model.params.segment[i]*model.params.b[i] for i in model.components)/1)]

function lb_volume(model::SAFTModel, z = SA[1.0]; phase = :unknown)
    seg = model.params.segment.values
    σᵢᵢ = model.params.sigma.diagvalues
    val = π/6*N_A*sum(z[i]*seg[i]*σᵢᵢ[i]^3 for i in 1:length(z))
    return val
end

function lb_volume(model::CubicModel,z = SA[1.0]; phase = :unknown)
    n = sum(z)
    invn = one(n)/n
    b = model.params.b.values
    b̄ = dot(z,Symmetric(b),z)*invn*invn
    return b̄
end

function lb_volume(model::CPAModel,z = SA[1.0]; phase = :unknown)
    n = sum(z)
    invn = one(n)/n
    b = model.params.b.values
    b̄ = dot(z,Symmetric(b),z)*invn*invn
    return b̄

end

function lb_volume(model::SAFTgammaMieModel, z = SA[1.0]; phase = :unknown)
    vk  = model.igroups
    seg = model.params.segment.values
    S   = model.params.shapefactor.values
    σᵢᵢ = model.params.sigma.diagvalues
    val = π/6*N_A*sum(z[i]*sum(vk[i][k]*seg[k]*S[k]*σᵢᵢ[k]^3 for k in @groups(i)) for i in @comps)
    return val
end

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



"""
    scale_sat_pure(model,T,z=SA[1.0])
returns the first guesses for the equilibrium volumes at a certain T



"""
function x0_sat_pure(model::EoSModel,T,z=SA[1.0])
    b = lb_volume(model,z)*one(T)
    B = second_virial_coefficient(model,T,z)
    x0v = -2*B + 2*b
    p = -0.25*R̄*T/B
    x0l = volume_compress(model,p,T,z)
    
    #=here we solve the saturation with aproximate 
    models of the EoS. on the gas side, we use
    a virial for logϕ, on the liquid side, we
    use an isothermal compressibility model based
    at the point P = P(-2B), vl(P)
   
    
    β = vt_isothermal_compressibility(model,x0l,T,z)
    fugv(P) = B*P/(R̄*T) #log(ϕv)
    P0 = p
    v0 = x0l
    fuglx(P) = -((v0*exp(β*P0)/R̄*T*β)*exp(-β*P) -log(P))
    fugl0 = fuglx(P0)
    fugl(P) = fuglx(P) - fugl0
    @show fugv(0.8P0)
    @show fugl(0.8P0)
    @show P0
    f0(P) = fugl(P) - fugv(P)
    @show Roots.find_zero(f0,0.9P0)
     =#
    x0  = [x0l,x0v]
    
    return log10.(x0)
end
#=
appendix: logϕ for the isothermal compressibility aprox
from the volume_compress code:

    Δ(P) =  -(P-p0)*β = -βP + βP0
    v(P) = v0*exp(Δ(P)) #volume_compress uses a convergence modification

    Z-1 = v/RT - 1

    definition of fugacity coefficient

    logϕ = ∫(Z-1)/P dp from Psat to Pmax
    Z(P) = P*v(P)/RT
    logϕ = ∫v/RT -1/P dP
    logϕ = ∫v0*exp(βP0)*exp(-βP)/RT -1/P dP
    logϕ = (v0*exp(βP0)/RT)∫exp(-βP) dp -  ∫1/P dp
    logϕ = -(v0*exp(βP0)/RTβ)*exp(-βP) -log(P)



=#

function scale_sat_pure(model::EoSModel,z=SA[1.0])
    p    = 1/p_scale(model,z)
    μ    = 1/R̄/T_scale(model,z)
    return p,μ
end


#=x0_crit_pure=#

# x0_crit_pure(model::LJSAFT) = [1.5, log10(π/6*model.params.segment[model.components[1]]*model.params.b[model.components[1]]/0.3)]

function x0_crit_pure(model::SAFTModel,z=SA[1.0])
    lb_v = lb_volume(model,z)
    [2, log10(lb_v/0.3)]
end

function x0_crit_pure(model::CubicModel,z=SA[1.0])
    lb_v = lb_volume(model,z)
    [1.0, log10(lb_v/0.3)]
end

function x0_crit_pure(model::CPAModel,z=SA[1.0])
    lb_v = lb_volume(model,z)
    [2.0, log10(lb_v/0.3)]
end

function x0_crit_pure(model::EoSModel,z=SA[1.0])
    lb_v = lb_volume(model,z)
    [1.5, log10(lb_v/0.3)]
end


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

T_crit_pure(model::EoSModel,z=SA[1.0]) = T_scale(model,z)

# T_crit_pure(model::Cubic) = model.params.a[model.components[1]]/model.params.b[model.components[1]]/8.314*8/27


#=
 temperature scaling factor,
on critical based EoS, is a function of critical temperature
on SAFT EoS, is a function of ϵ
=#
function T_scale(model::SAFTModel,z=SA[1.0])
    ϵ = model.params.epsilon.diagvalues
    return prod(ϵ)^(1/length(ϵ))
end

#dont use αa, just a, to avoid temperature dependence
function T_scale(model::CubicModel,z=SA[1.0])
    n = sum(z)
    invn2 = one(n)/(n*n)
    Ωa,Ωb = ab_consts(model)
    _a = model.params.a.values
    _b = model.params.b.values
    a = dot(z, Symmetric(_a), z)*invn2/Ωa
    b = dot(z, Symmetric(_b), z)*invn2/Ωb
    return a/b/R̄
end

function T_scale(model::CPAModel,z=SA[1.0])
    n = sum(z)
    invn2 = one(n)/(n*n)
    Ωa,Ωb = ab_consts(model)
    _a = model.params.a.values
    _b = model.params.b.values
    a = dot(z, Symmetric(_a), z)*invn2/Ωa
    b = dot(z, Symmetric(_b), z)*invn2/Ωb
    return a/b/R̄

end

#=
pressure scaling factor
on critical eos, a function of critical pressure
on SAFT, a function of
=#
function p_scale(model::SAFTModel,z=SA[1.0])
    ϵ = model.params.epsilon.diagvalues
    σᵢᵢ = model.params.sigma.diagvalues
    val =  sum(z[i]*σᵢᵢ[i]^3/ϵ[i] for i in 1:length(z))*N_A/R̄
    return 1/val
end

function p_scale(model::SAFTgammaMieModel,z=SA[1.0])
    vk  = model.igroups
    seg = model.params.segment.values
    S   = model.params.shapefactor.values
    σ   = model.params.sigma.values
    m̄  = sum(vk[1][k]*seg[k]*S[k] for k in @groups(model.icomponents[1]))

    σ̄3 = sum(vk[1][k]*vk[1][l]*
                     seg[k]*seg[l]*
                     S[k]*S[l]*
                     σ[k,l]^3 for k in @groups(model.icomponents[1]) for l in @groups(model.icomponents[1]))/m̄^2
    ϵ̄ = T_scale(model,z)
    val    = σ̄3*N_A/R̄/ϵ̄
    return 1/val
end

function p_scale(model::CubicModel,z=SA[1.0])
    n = sum(z)
    invn2 = (1/n)^2
    Ωa,Ωb = ab_consts(model)
    _a = model.params.a.values
    _b = model.params.b.values
    a = invn2*dot(z, Symmetric(_a), z)/Ωa
    b = invn2*dot(z, Symmetric(_b), z)/Ωb
    return a/ (b^2) # Pc mean
end

function p_scale(model::CPAModel,z=SA[1.0])
    n = sum(z)
    invn2 = (1/n)^2
    Ωa,Ωb = ab_consts(model)
    _a = model.params.a.values
    _b = model.params.b.values
    a = invn2*dot(z, Symmetric(_a), z)/Ωa
    b = invn2*dot(z, Symmetric(_b), z)/Ωb
    return a/ (b^2) # Pc mean
end



#=
the following methods are fallbacks,
that require just the definition of T_scale,p_scale and lb_volume
respectively. if possible, each eos should define those
=#
function T_scales(model,z)
    n = length(z)
    x = zeros(n)
    res = zeros(n)
    for i = 1:n
        x[i] = 1.0
        res[i] = T_scale(model,x)
        x[i] = 0.0
    end
    return res
end


function lb_volumes(model,z)
    n = length(z)
    x = zeros(n)
    res = zeros(n)
    for i = 1:n
        x[i] = 1.0
        res[i] = lb_volume(model,x)
        x[i] = 0.0
    end
    return res
end

function p_scales(model,z)
    n = length(z)
    x = zeros(n)
    res = zeros(n)
    for i = 1:n
        x[i] = 1.0
        res[i] = p_scale(model,x)
        x[i] = 0.0
    end
end
