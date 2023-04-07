"""
    __gas_model(model::EoSModel)

internal function.
provides the model used to calculate gas properties.
normally, this is the identity, but `CompositeModel` has a gas model by itself.
"""
__gas_model(model::EoSModel) = model

include("bubble_point.jl")
include("dew_point.jl")
include("LLE_point.jl")




function PTFlashWrapper(model::ActivityModel,T::Number) 
    pures = model.puremodel.pure
    sats = saturation_pressure.(pures,T)
    vv_pure = last.(sats)
    RT = R̄*T
    p_pure = first.(sats)
    μpure = only.(VT_chemical_potential_res.(__gas_model.(pures),vv_pure,T))
    ϕpure = exp.(μpure ./ RT .- log.(p_pure .* vv_pure ./ RT))
    return PTFlashWrapper(model.components,model,sats,ϕpure)
end

__tpflash_cache_model(model::ActivityModel,p,T,z) = PTFlashWrapper(model,T)

function update_K!(lnK,wrapper::PTFlashWrapper{<:ActivityModel},p,T,x,y,volx,voly,phasex,phasey,β = nothing,inx = FillArrays.Fill(true,length(x)),iny = inx)
    model = wrapper.model
    pures = wrapper.model.puremodel.pure
    sats = wrapper.sat
    n = length(model)
    #crits = wrapper.crit
    fug = wrapper.fug
    RT = R̄*T
    γx = activity_coefficient(model, p, T, x)
    volx = volume(model.puremodel.model, p, T, x, phase = phasex, vol0 = volx)
    _0 = zero(eltype(lnK))

    if β === nothing
        _0 = zero(eltype(lnK))
        gibbs = _0/_0
    else
        gibbs = _0
        for i in eachindex(x)
            if inx[i]
                g_E_x = x[i]*RT*log(γx[i])
                g_ideal_x = x[i]*RT*log(x[i])
                g_pure_x = x[i]*VT_gibbs_free_energy(pures[i],wrapper.sat[i][2],T)
                gibbs += (g_E_x + g_ideal_x + g_pure_x)*(1-β)/RT
            end
        end
    end
    
    if is_vapour(phasey)
        lnϕy, voly = lnϕ(model, p, T, y; phase=phasey, vol0=voly)
        for i in eachindex(lnK)
            if iny[i]
                ϕli = fug[i]
                p_i = sats[i][1]
                lnK[i] = log(γx[i]*p_i*ϕli/p) - lnϕy[i] + volx*(p - p_i)/RT
                gibbs += β*y[i]*log(y[i] + lnϕy[i])
            end
        end
    else
        γy = activity_coefficient(model, p, T, y)
        lnK .= log.(γx./γy)
        voly = volume(model.puremodel.model, p, T, y, phase = phasey, vol0 = voly)
        if β !== nothing
            for i in eachindex(y)
                if iny[i]
                    g_E_y = y[i]*RT*log(γy[i])
                    g_ideal_y = y[i]*RT*(log(y[i]))
                    g_pure_y = y[i]*VT_gibbs_free_energy(pures[i],wrapper.sat[i][2],T)
                    gibbs += (g_E_y + g_ideal_y + g_pure_y)*β/RT
                end
            end
        end
    end
    
    return lnK,volx,voly,gibbs
end

function __tpflash_gibbs_reduced(wrapper::PTFlashWrapper{<:ActivityModel},p,T,x,y,β,eq)
    pures = wrapper.model.puremodel.pure
    model = wrapper.model
    γx = activity_coefficient(model, p, T, x)
    RT = R̄*T
    n = length(model)
    g_E_x = sum(x[i]*RT*log(γx[i]) for i ∈ 1:n)
    g_ideal_x = sum(x[i]*RT*(log(x[i])) for i ∈ 1:n)
    g_pure_x = sum(x[i]*VT_gibbs_free_energy(pures[i],wrapper.sat[i][2],T) for i ∈ 1:n)    
    gibbs = (g_E_x + g_ideal_x + g_pure_x)*(1-β)/RT
    if is_vle(eq)
        gibbs += gibbs_free_energy(model.puremodel.model,p,T,y,phase =:v)*β/R̄/T
    else #lle
        γy = activity_coefficient(model, p, T, y)
        g_E_y = sum(y[i]*RT*log(γy[i]) for i ∈ 1:n)
        g_ideal_y = sum(y[i]*R̄*T*(log(y[i])) for i ∈ 1:n)
        g_pure_y = sum(y[i]*VT_gibbs_free_energy(pures[i],wrapper.sat[i][2],T) for i ∈ 1:n)    
        gibbs += (g_E_y + g_ideal_y + g_pure_y)*β/RT
    end
    return gibbs
    #(gibbs_free_energy(model,p,T,x)*(1-β)+gibbs_free_energy(model,p,T,y)*β)/R̄/T
end

function dgibbs_obj!(model::PTFlashWrapper{<:ActivityModel}, p, T, z, phasex, phasey,
    nx, ny, vcache, ny_var = nothing, in_equilibria = FillArrays.Fill(true,length(z)), non_inx = in_equilibria, non_iny = in_equilibria;
    F=nothing, G=nothing, H=nothing)

    # Objetive Function to minimize the Gibbs Free Energy
    # It computes the Gibbs free energy, its gradient and its hessian
    iv = 0
    for i in eachindex(z)
        if in_equilibria[i]
            iv += 1
            nyi = ny_var[iv]
            ny[i] = nyi
            nx[i] =z[i] - nyi
        end
    end    # nx = z .- ny

    nxsum = sum(nx)
    nysum = sum(ny)
    x = nx ./ nxsum
    y = ny ./ nysum

    # Volumes are set from local cache to reuse their values for following
    # Iterations
    volx,voly = vcache[]
    all_equilibria = all(in_equilibria)
    if H !== nothing
        # Computing Gibbs Energy Hessian
        lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, volx = ∂lnϕ∂n∂P(model, p, T, x; phase=phasex, vol0=volx)
        lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, voly = ∂lnϕ∂n∂P(model, p, T, y; phase=phasey, vol0=voly)

        if !all_equilibria
            ∂ϕx = ∂lnϕ∂nx[in_equilibria, in_equilibria]
            ∂ϕy = ∂lnϕ∂ny[in_equilibria, in_equilibria]
        else
            #skip a copy if possible
            ∂ϕx,∂ϕy = ∂lnϕ∂nx,∂lnϕ∂ny
        end
            ∂ϕx .-= 1
            ∂ϕy .-= 1
            ∂ϕx ./= nxsum
            ∂ϕy ./= nysum
        for (i,idiag) in pairs(diagind(∂ϕy))
            ∂ϕx[idiag] += 1/nx[i]
            ∂ϕy[idiag] += 1/ny[i]
        end

        #∂ϕx = eye./nx .- 1/nxsum .+ ∂lnϕ∂nx/nxsum
        #∂ϕy = eye./ny .- 1/nysum .+ ∂lnϕ∂ny/nysum
        H .= ∂ϕx .+ ∂ϕy
    else
        lnϕx, volx = lnϕ(model, p, T, x; phase=phasex, vol0=volx)
        lnϕy, voly = lnϕ(model, p, T, y; phase=phasey, vol0=voly)
    end
    #volumes are stored in the local cache
    vcache[] = (volx,voly)

    ϕx = log.(x) .+ lnϕx
    ϕy = log.(y) .+ lnϕy

    # to avoid NaN in Gibbs energy
    for i in eachindex(z)
        non_iny[i] && (ϕy[i] = 0.)
        non_inx[i] && (ϕx[i] = 0.)
    end

    if G !== nothing
        # Computing Gibbs Energy gradient
        if !all_equilibria
            G .= (ϕy .- ϕx)[in_equilibria]
        else
            G .= ϕy .- ϕx
        end
    end

    if F !== nothing
        # Computing Gibbs Energy
        FO = dot(ny,ϕy) + dot(nx,ϕx)
        return FO
    end
end
