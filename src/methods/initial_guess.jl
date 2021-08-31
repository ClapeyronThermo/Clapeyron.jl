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

"""
    x0_volume_liquid(model,T,z)

Returns an initial guess to the liquid volume, dependent on temperature and composition. by default is 1.25 times the [lower bound volume](@ref lb_volume).
"""
function x0_volume_liquid(model,T,z)
    v_lb = lb_volume(model,z)
    return v_lb*1.25
end

function x0_volume_liquid(model::SAFTVRMieModel,T,z)
    v_lb = lb_volume(model,z)
    return v_lb*1.5
end

function x0_volume_liquid(model::SAFTgammaMieModel,T,z)
    v_lb = lb_volume(model,z)
    return v_lb*1.5
end
"""
    x0_volume_gas(model,p,T,z)

Returns an initial guess to the gas volume, depending of pressure, temperature and composition. by default uses a [virial aproximation](@ref volume_virial)
"""
function x0_volume_gas(model,p,T,z)
    return volume_virial(model,p,T,z)
end

function x0_volume_sc(model,p,T,z)
    v_sc = lb_volume(model,z)
    return v_sc*2
end

"""
    x0_volume(model::EoSModel,p,T,z; phase = :unknown)

Returns an initial guess of the volume at a pressure, temperature, composition and suggested phase.

If the suggested phase is `:unkwown` or `:liquid`, calls [`x0_volume_liquid`](@ref).

If the suggested phase is `:gas`, calls [`x0_volume_gas`](@ref).

"""
function x0_volume(model::EoSModel,p,T,z; phase = :unknown)
    phase = Symbol(phase)
    if phase === :unknown || is_liquid(phase)
        return x0_volume_liquid(model,T,z)
    elseif is_vapour(phase)
        return x0_volume_gas(model,p,T,z)
    elseif is_supercritical(phase)
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



"""
    lb_volume(model::EoSModel,z=SA[1.0])

Returns the lower bound volume. 

It has different meanings depending on the Equation of State, but symbolizes the minimum allowable volume at a certain composition:

- SAFT EoS: the packing volume
- Cubic EoS, covolume (b) parameter

On empiric equations of state, the value is chosen to match the volume of the conditions at maximum pressure and minimum temperature
, but the equation itself normally can be evaluated at lower volumes.

On SAFT and Cubic EoS, volumes lower than `lb_volume` will likely error.

The lower bound volume is used for guesses of liquid volumes at a certain pressure, saturated liquid volumes and critical volumes.
"""
function lb_volume end


function lb_volume(model::SAFTModel, z = SA[1.0])
    seg = model.params.segment.values
    σᵢᵢ = model.params.sigma.diagvalues
    val = π/6*N_A*sum(z[i]*seg[i]*σᵢᵢ[i]^3 for i in 1:length(z))
    return val
end

function lb_volume(model::CubicModel,z = SA[1.0])
    n = sum(z)
    V = 1e-5
    T = 0.
    invn = one(n)/n
    b = model.params.b.values
    c = @f(translation,model.translation)
    b̄ = dot(z,Symmetric(b),z)*invn*invn
    c̄ = dot(z,c)*invn
    return b̄-c̄
end

function lb_volume(model::CPAModel,z = SA[1.0])
    n = sum(z)
    invn = one(n)/n
    b = model.params.b.values
    b̄ = dot(z,Symmetric(b),z)*invn*invn
    return b̄

end

function lb_volume(model::SAFTgammaMieModel, z = SA[1.0])
    vk  = model.groups.n_flattenedgroups
    seg = model.params.segment.values
    S   = model.params.shapefactor.values
    σᵢᵢ = model.params.sigma.diagvalues
    val = π/6*N_A*sum(z[i]*sum(vk[i][k]*seg[k]*S[k]*σᵢᵢ[k]^3 for k in @groups(i)) for i in @comps)
    return val
end

function lb_volume(model::BACKSAFTModel, z = SA[1.0])
    seg = model.params.segment.values
    σᵢᵢ = model.params.sigma.diagvalues
    α   = model.params.alpha.values
    val = π/6*N_A*sum(z[i]*α[i]*seg[i]*σᵢᵢ[i]^3 for i in 1:length(z))
    return val
end

function lb_volume(model::LJSAFTModel, z = SA[1.0])
    seg = model.params.segment.values
    b = model.params.b.diagvalues
    val = π/6*sum(z[i]*seg[i]*b[i] for i in 1:length(z))
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
    x0_sat_pure(model::EoSModel,T,z=SA[1.0])

Returns a 2-element vector corresponding to `[log10(Vₗ),log10(Vᵥ)]`, where Vₗ and Vᵥ are the liquid and vapor initial guesses. 
Used in [`sat_pure`](@ref).
"""
function x0_sat_pure(model,T,z=SA[1.0])
    #=theory as follows
    #given T = Teos:
    #calculate B(Teos)
    PB = Peos = -2B
    
    veos = volume_compress(model,PB,T)
    with a (P,V,T,B) pair, change to
    (P,V,a,b)
=#
    B = second_virial_coefficient(model,T,SA[1.0])
    lb_v = lb_volume(model,SA[1.0])*one(T)
    _0 = zero(B)
    #virial volume below lower bound volume.
    #that means that we are way over the critical point
    if -2B < lb_v 
        _nan = _0/_0
        return [_nan,_nan]
    end
    p = -0.25*R̄*T/B
    vl = x0_volume(model,p,T,z,phase=:l)
    
    vl = volume_compress(model,p,T,z,V0=vl)
    #=the basis is that p = RT/v-b - a/v2
    we have a (p,v,T) pair
    and B = 2nd virial coefficient = b-a/RT
    with that, we solve for a and b 
    as a and b are vdW, Pc and Tc can be calculated
    with Tc and Pc, we use [1] to calculate vl0 and vv0
    with Tc, we can also know in what regime we are.
    in near critical pressures, we use directly vv0 = -2B
    and vl0 = 4*lb_v

    [1]
    DOI: 10.1007/s10910-007-9272-4
    Journal of Mathematical Chemistry, Vol. 43, No. 4, May 2008 (© 2007)
    The van der Waals equation: analytical and approximate solutions
    =#
    γ = p*vl*vl/(R̄*T)
    _c = vl*vl + B*vl - γ*vl 
    _b = γ - B - vl
    Δ = _b*_b - 4*_c
    if isnan(vl) | (Δ < 0)
        #old strategy
        x0l = 4*lb_v
        x0v = -2*B + 2*lb_v
        return [log10(x0l),log10(x0v)]   
    end
    Δsqrt = sqrt(Δ)
    b1 = 0.5*(-_b + Δsqrt)
    b2 = 0.5*(-_b - Δsqrt)
    if b1 < 0
        b = b2
    elseif b1 > vl
        b = b2
    else
        b = b1
    end
    a = -R̄*T*(B-b)
    Vc = 3*b
    Ωa =  27/64
    Ωb =  1/8
    ar = a/Ωa
    br = b/Ωb
    Tc = ar/br/R̄
    Pc = ar/(br*br)
    Tr = T/Tc
    #Tr(vdW approx) > Tr(model)
    if Tr >= 1
        _nan = _0/_0
        return [_nan,_nan]
    end
    # if b1/b2 < -0.95, then T is near Tc.
    #if b<lb_v then we are in trouble 
    #critical regime or something horribly wrong happened
    if (b1/b2 < -0.95) | (b<lb_v) | (Tr>0.99)
        x0l = 4*lb_v
        x0v = -2*B #gas volume as high as possible
        return [log10(x0l),log10(x0v)]   
    end
   
    Trm1 = 1.0-Tr
    Trmid = sqrt(Trm1)
    if Tr >= 0.7 
        c_l = 1.0+2.0*Trmid + 0.4*Trm1 - 0.52*Trmid*Trm1 +0.115*Trm1*Trm1 #Eq. 29
    else
        c_l = 1.5*(1+sqrt(1-(32/27)*Tr)) #Eq. 32
    end

    if Tr >= 0.46
        #Eq. 30, valid in 0.46 < Tr < 1 
        c_v = 1.0-2.0*Trmid + 0.4*Trm1 + 0.52*Trmid*Trm1 +0.207*Trm1*Trm1   
    elseif Tr <= 0.33
        #Eq. 33, valid in 0 < Tr < 0.33
        c_v = (3*c_l/(ℯ*(3-c_l)))*exp(-(1/(1-c_l/3)))

    else
        #Eq. 31 valid in 0.25 < Tr < 1
        mean_c = 1.0 + 0.4*Trm1 + 0.161*Trm1*Trm1     
        c_v = 2*mean_c - c_l

    end
    #volumes predicted by vdW
    Vl0 = (1/c_l)*Vc
    Vv0 = (1/c_v)*Vc 
    x0l = min(Vl0,vl)
    x0v = min(1e4*one(Vv0),Vv0) #cutoff volume
    return [log10(x0l),log10(x0v)] 
end

function scale_sat_pure(model::EoSModel,z=SA[1.0])
    p    = 1/p_scale(model,z)
    μ    = 1/R̄/T_scale(model,z)
    return p,μ
end


#=x0_crit_pure=#

# x0_crit_pure(model::LJSAFT) = [1.5, log10(π/6*model.params.segment[model.components[1]]*model.params.b[model.components[1]]/0.3)]

"""
    x0_crit_pure(model::SAFTModel,z=SA[1.0])

Returns a 2-element vector corresponding to
    `[k,log10(Vc0)]`, where `k` is `Tc0/T_scale(model,z)`
"""
function x0_crit_pure end

function x0_crit_pure(model::SAFTModel,z=SA[1.0])
    lb_v = lb_volume(model,z)
    [2.0, log10(lb_v/0.3)]
end

function x0_crit_pure(model::sCKSAFTModel,z=SA[1.0])
    lb_v = lb_volume(model,z)
    res = [5.0, log10(lb_v/0.3)]
    #@show res
    return res
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


# T_crit_pure(model::Cubic) = model.params.a[model.components[1]]/model.params.b[model.components[1]]/8.314*8/27


#=
 temperature scaling factor,
on critical based EoS, is a function of critical temperature
on SAFT EoS, is a function of ϵ
=#


"""
    T_scale(model::EoS,z=SA[1.0])

Represents a temperature scaling factor. 

On any EoS based on Critical parameters (Cubic or Empiric EoS), the temperature scaling factor is chosen to be the critical temperature.

On SAFT or other molecular EoS, the temperature scaling factor is chosen to be a function of the potential depth ϵ.

Used as scaling factors in [`sat_pure`](@ref) and as input for solving [`crit_pure`](@ref)
"""
function T_scale end

function T_scale(model::SAFTModel,z=SA[1.0])
    ϵ = model.params.epsilon.diagvalues
    return prod(ϵ)^(1/length(ϵ))
end

function T_scale(model::LJSAFTModel,z=SA[1.0])
    T̃ = model.params.T_tilde.diagvalues
    return prod(T̃)^(1/length(T̃))
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

"""
    p_scale(model::SAFTModel,z=SA[1.0])

Represents a pressure scaling factor

On any EoS based on Critical parameters (Cubic or  
Empiric EoS), the pressure scaling factor is    
chosen to be a function of the critical pressure.

On SAFT or other molecular EoS, the temperature    
scaling factor is chosen to a function of ∑(zᵢ*ϵᵢ*(σᵢᵢ)³)    

Used as scaling factors in [`sat_pure`](@ref) and as input for solving [`crit_pure`](@ref)

"""
function p_scale(model::SAFTModel,z=SA[1.0])
    ϵ = model.params.epsilon.diagvalues
    σᵢᵢ = model.params.sigma.diagvalues
    val =  sum(z[i]*σᵢᵢ[i]^3/ϵ[i] for i in 1:length(z))*N_A/R̄
    return 1/val
end

function p_scale(model::LJSAFTModel,z=SA[1.0])
    T̃ = model.params.T_tilde.diagvalues
    b = model.params.b.diagvalues
    val =  sum(z[i]*b[i]/T̃[i] for i in 1:length(z))/R̄
    return 1/val
end

function p_scale(model::SAFTgammaMieModel,z=SA[1.0])
    vk  = model.groups.n_flattenedgroups
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


# function lb_volumes(model,z)
#     n = length(z)
#     x = zeros(n)
#     res = zeros(n)
#     for i = 1:n
#         x[i] = 1.0
#         res[i] = lb_volume(model,x)
#         x[i] = 0.0
#     end
#     return res
# end

# function p_scales(model,z)
#     n = length(z)
#     x = zeros(n)
#     res = zeros(n)
#     for i = 1:n
#         x[i] = 1.0
#         res[i] = p_scale(model,x)
#         x[i] = 0.0
#     end
# end