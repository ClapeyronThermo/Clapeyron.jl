#special support for activity models with the michelsen tp flash method.

function PTFlashWrapper(model::ActivityModel,T::Number) 
    pures = model.puremodel.pure
    sats = saturation_pressure.(pures,T)
    vv_pure = last.(sats)
    RT = R̄*T
    p_pure = first.(sats)
    μpure = only.(VT_chemical_potential_res.(__gas_model.(pures),vv_pure,T))
    ϕpure = exp.(μpure ./ RT .- log.(p_pure .* vv_pure ./ RT))
    g_pure = [VT_gibbs_free_energy(__gas_model(pures[i]),sats[i][2],T) for i in 1:length(model)]
    return PTFlashWrapper(model.components,model,sats,ϕpure,g_pure)
end

__tpflash_cache_model(model::ActivityModel,p,T,z) = PTFlashWrapper(model,T)

function update_K!(lnK,wrapper::PTFlashWrapper{<:ActivityModel},p,T,x,y,volx,voly,phasex,phasey,β = nothing,inx = FillArrays.Fill(true,length(x)),iny = inx)
    model = wrapper.model
    sats = wrapper.sat
    g_pures = wrapper.μ
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
                g_pure_x = x[i]*g_pures[i]
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
                    g_pure_y = y[i]*g_pures[i]
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
    g_pures = wrapper.μ

    γx = activity_coefficient(model, p, T, x)
    RT = R̄*T
    n = length(model)
    g_E_x = sum(x[i]*RT*log(γx[i]) for i ∈ 1:n)
    g_ideal_x = sum(x[i]*RT*(log(x[i])) for i ∈ 1:n)
    g_pure_x = dot(x,g_pures)
    gibbs = (g_E_x + g_ideal_x + g_pure_x)*(1-β)/RT
    if is_vle(eq)
        gibbs += gibbs_free_energy(model.puremodel.model,p,T,y,phase =:v)*β/R̄/T
    else #lle
        γy = activity_coefficient(model, p, T, y)
        g_E_y = sum(y[i]*RT*log(γy[i]) for i ∈ 1:n)
        g_ideal_y = sum(y[i]*R̄*T*(log(y[i])) for i ∈ 1:n)
        g_pure_y = dot(y,g_pures)
        gibbs += (g_E_y + g_ideal_y + g_pure_y)*β/RT
    end
    return gibbs
    #(gibbs_free_energy(model,p,T,x)*(1-β)+gibbs_free_energy(model,p,T,y)*β)/R̄/T
end

#TODO: derive expressions for this

function dgibbs_obj!(model::PTFlashWrapper{<:ActivityModel}, p, T, z, phasex, phasey,
    nx, ny, vcache, ny_var = nothing, in_equilibria = FillArrays.Fill(true,length(z)), non_inx = in_equilibria, non_iny = in_equilibria;
    F=nothing, G=nothing, H=nothing)
    throw(error("Activity models don't support gibbs energy optimization in MichelsenTPFlash."))

end