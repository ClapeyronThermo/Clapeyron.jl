"""
    GammaPhi{γ,Φ} <: RestrictedEquilibriaModel

Wrapper struct to signal that a `CompositeModel` uses an activity model in conjunction with a fluid.
"""
struct GammaPhi{γ,Φ} <: RestrictedEquilibriaModel
    components::Vector{String}
    activity::γ
    fluid::EoSVectorParam{Φ}
end

reference_state(model::GammaPhi) = reference_state(model.fluid)

function Base.show(io::IO,mime::MIME"text/plain",model::GammaPhi)
    print(io,"γ-ϕ Model")
    length(model) == 1 && print(io, " with 1 component:")
    length(model) > 1 && print(io, " with ", length(model), " components:")
    println(io)
    show_pairs(io,model.components)
    act = model.activity
    if hasfield(typeof(act),:puremodel)
        print(io,'\n',"Activity Model: ", parameterless_type(act))
    else
        print(io,'\n',"Activity Model: ",typeof(act))
    end
    print(io,'\n',"Fluid Model: ",typeof(model.fluid.model))
    show_reference_state(io,model;space = true)
end

fluid_model(model::GammaPhi) = model.fluid.model

function activity_coefficient(model::GammaPhi,p,T,z=SA[1.];
                            μ_ref = nothing,
                            reference = :pure,
                            phase=:unknown,
                            threaded=true,
                            vol0=nothing)

    return activity_coefficient(model.activity,p,T,z;μ_ref,reference,phase,threaded,vol0)
end

function excess_gibbs_free_energy(model::GammaPhi,p,T,z)
    return excess_gibbs_free_energy(model.activity,p,T,z)
end

reference_chemical_potential_type(model::GammaPhi) = reference_chemical_potential_type(model.activity)

function volume_impl(model::GammaPhi,p,T,z,phase,threaded,vol0)
    return volume_impl(model.fluid.model,p,T,z,phase,threaded,vol0)
end

molecular_weight(model::GammaPhi,z) = molecular_weight(model.fluid.model,z)
saturation_model(model::GammaPhi) = saturation_model(model.fluid)
idealmodel(model::GammaPhi) = idealmodel(model.fluid.model)

function init_preferred_method(method::typeof(tp_flash),model::GammaPhi,kwargs)
    RRTPFlash(;kwargs...)
end

# Error handling for Activity models that don't provide saturation properties, in the context of VLE.
function ActivitySaturationError(model,method)
    throw(ArgumentError("$method requires $model to be used in conjuction with another EoS model that supports saturation properties. If you are using an Activity Model as a raw input, use `CompositeModel(components, liquid = activity_model, fluid = fluid_model)` instead."))
end

struct ActivityEval{γ} <: IdealModel
    gammaphi::γ
end

function gibbs_solvation(model::GammaPhi,T)
    z = [1.0,1e-30]
    p,v_l,v_v = saturation_pressure(model.fluid[1],T)
    p2,v_l2,v_v2 = saturation_pressure(model.fluid[2],T)
    γ = activity_coefficient(model,p,T,z)
    K = v_v/v_l*γ[2]*p2/p
    return -R̄*T*log(K)
end


#=
some equations require more delicate handling. the general case only works for extensive properties.
=#
function PT_property_gammaphi(model::GammaPhi,p,T,z,f::F,USEP) where F
    ∑z = sum(z)
    z1 = SA[∑z]
    res = zero(Base.promote_eltype(model,p,T,z))
    Vl = zero(res)
    use_p = USEP === Val{true}()
    for i in 1:length(model)
        zi = z[i]
        modeli = model.fluid.pure[i]
        #we suppose this is liquid phase
        Vi = volume(modeli, p, T, z1; phase = :liquid,threaded = false)
        Vl += Vi*zi/∑z
        if use_p
            res += f(modeli,Vi,T,z1,p)*zi/∑z
        else
            res += f(modeli,Vi,T,z1)*zi/∑z
        end
    end
    if !iszero(p)
        if use_p
            res += f(model.activity,0.0,T,z,p)
        else
            res += f(model.activity,0.0,T,z)
        end
        return res
    end
end

function gammaphi_f_hess(model,p,T,z)
    ∑z = sum(z)
    z1 = SA[∑z]
    _0 = zero(Base.promote_eltype(model,p,T,z))
    _1 = oneunit(_0)
    Vl = zero(_0)
    ∂²A∂T² = ∂²f∂T²(model.activity,_0,T,z)
    ∂²A = SMatrix{2,2}((_0,_0,_0,_0 + ∂²A∂T²))
    for i in 1:length(model)
        zi = z[i]
        modeli = model.fluid.pure[i]
        #eos(model) = sum(xi*eos(model,Vi,T,1.0))

        Vi = volume(modeli, p, T, z1; phase = :liquid, threaded = false)
        Vl += Vi*zi/∑z
        Vii = Vi*zi
        ∂²Aᵢ = f_hess(modeli,Vi,T,z1) .* zi ./ ∑z
        ∂²A = ∂²A .+ ∂²Aᵢ
    end
    return ∂²A,Vl
end

function PT_property_gammaphi(model::GammaPhi,p,T,z,::typeof(VT_isobaric_heat_capacity),USEP)
    d²A,V = gammaphi_f_hess(model,p,T,z)
    ∂²A∂V∂T = d²A[1,2]
    ∂²A∂V² = d²A[1,1]
    ∂²A∂T² = d²A[2,2]
    return -T*(∂²A∂T² - ∂²A∂V∂T^2/∂²A∂V²)
end

function PT_property_gammaphi(model::GammaPhi,p,T,z,::typeof(VT_adiabatic_index),USEP)
    d²A,V = gammaphi_f_hess(model,p,T,z)
    ∂²A∂V∂T = d²A[1,2]
    ∂²A∂V² = d²A[1,1]
    ∂²A∂T² = d²A[2,2]
    return 1 - ∂²A∂V∂T*∂²A∂V∂T/(∂²A∂V²*∂²A∂T²)
end

function PT_property_gammaphi(model::GammaPhi,p,T,z,::typeof(VT_isothermal_compressibility),USEP)
    d²A,V = gammaphi_f_hess(model,p,T,z)
    return -1/V/d²A[1,1]
end

function PT_property_gammaphi(model::GammaPhi,p,T,z,::typeof(VT_isentropic_compressibility),USEP)
    d²A,V = gammaphi_f_hess(model,p,T,z)
    ∂²A∂V∂T = d²A[1,2]
    ∂²A∂V² = d²A[1,1]
    ∂²A∂T² = d²A[2,2]
    return 1/V/(∂²A∂V²-∂²A∂V∂T^2/∂²A∂T²)
end



function PT_property_gammaphi(model::GammaPhi,p,T,z,::typeof(VT_speed_of_sound),USEP)
    Mr = molecular_weight(model,z)
    d²A,V = gammaphi_f_hess(model,p,T,z)
    ∂²A∂V∂T = d²A[1,2]
    ∂²A∂V² = d²A[1,1]
    ∂²A∂T² = d²A[2,2]
    return V*sqrt((∂²A∂V²-∂²A∂V∂T^2/∂²A∂T²)/Mr)
end

function PT_property_gammaphi(model::GammaPhi,p,T,z,::typeof(VT_isobaric_expansivity),USEP)
    d²A,V = gammaphi_f_hess(model,p,T,z)
    ∂²A∂V∂T = d²A[1,2]
    ∂²A∂V² = d²A[1,1]
    return -∂²A∂V∂T/(V*∂²A∂V²)
end


function PT_property_gammaphi(model::GammaPhi,p,T,z,::typeof(VT_joule_thomson_coefficient),USEP)
    d²A,V = gammaphi_f_hess(model,p,T,z)
    ∂²A∂V∂T = d²A[1,2]
    ∂²A∂V² = d²A[1,1]
    ∂²A∂T² = d²A[2,2]
    return -(∂²A∂V∂T - ∂²A∂V²*((T*∂²A∂T² + V*∂²A∂V∂T) / (T*∂²A∂V∂T + V*∂²A∂V²)))^-1
end

function PT_property(model::GammaPhi,p,T,z,phase,threaded,vol0,f::F,USEP::Val{UseP}) where {F,UseP}
    if phase == :stable
        #we dont support this in GammaPhi.
        throw(error("automatic phase detection not implemented for $(typeof(model))"))
    end

    #shortcut for one-component models:
    if length(model) == 1
        res_v = PT_property(model.fluid,p,T,z,phase,threaded,vol0,f,USEP)
    end
    #=
    Calculate PT properties as the following:
    prop = ∑xi*prop(pure[i],p,T) + excess_prop(activity,T,x) + prop(reference_state,T)

    Vapour properties are calculated with the fluid model
    =#
    if is_vapour(phase)
        res = PT_property(model.fluid,p,T,z,phase,threaded,vol0,f,USEP)
        return res
    elseif is_liquid(phase) || is_unknown(phase)
       return PT_property_gammaphi(model,p,T,z,f,USEP)
    else
        throw(error("invalid phase specifier: $phase"))
    end
end

function __calculate_reference_state_consts(model::GammaPhi,v,T,p,z,H0,S0,phase)
    ∑z = sum(z)
    S00 = entropy(model,p,T,z,phase = phase)
    a1 = (S00 - S0)#/∑z
    H00 = enthalpy(model,p,T,z,phase = phase)
    a0 = (-H00 + H0)#/∑z
    return a0,a1
end

function PTFlashWrapper(model::GammaPhi,p,T::Number,equilibrium::Symbol)
    fluidmodel = model.fluid
    #check that we can actually solve the equilibria
    if fluidmodel isa IdealModel && !is_lle(equilibrium)
        ActivitySaturationError(model.activity,tp_flash)
    end
    pures = fluidmodel.pure
    RT = R̄*T
    if fluidmodel.model isa IdealModel
        vv = RT/p
        nan = zero(vv)/zero(vv)
        sats = fill((nan,nan,vv),length(model))
        ϕpure = fill(one(vv),length(model))
        g_pure = [VT_gibbs_free_energy(gas_model(pures[i]),vv,T) for i in 1:length(model)]
        return PTFlashWrapper(model.components,model,sats,ϕpure,g_pure,equilibrium)
    else
        sats = saturation_pressure.(pures,T)
        vv_pure = last.(sats)
        p_pure = first.(sats)
        μpure = only.(VT_chemical_potential_res.(gas_model.(pures),vv_pure,T))
        ϕpure = exp.(μpure ./ RT .- log.(p_pure .* vv_pure ./ RT))
        g_pure = [VT_gibbs_free_energy(gas_model(pures[i]),vv_pure[i],T) for i in 1:length(model)]
        return PTFlashWrapper(model.components,model,sats,ϕpure,g_pure,equilibrium)
    end
end

__tpflash_cache_model(model::GammaPhi,p,T,z,equilibrium) = PTFlashWrapper(model,p,T,equilibrium)

function update_K!(lnK,wrapper::PTFlashWrapper{<:GammaPhi},p,T,x,y,z,β,vols,phases,non_inw,cache = nothing)
    volx,voly = vols
    phasex,phasey = phases
    non_inx,non_iny = non_inw
    model = wrapper.model
    fluidmodel = model.fluid.model
    sats = wrapper.sat
    g_pures = wrapper.μ
    n = length(model)
    #crits = wrapper.crit
    fug = wrapper.fug
    RT = R̄*T
    γx = activity_coefficient(model, p, T, x)
    volx = volume(fluidmodel, p, T, x, phase = phasex, vol0 = volx)
    _0 = zero(eltype(lnK))
    is_ideal = fluidmodel isa IdealModel
    if β === nothing
        _0 = zero(eltype(lnK))
        gibbs = _0/_0
    else
        gibbs = _0
        for i in eachindex(x)
            if !non_inx[i]
                g_E_x = x[i]*RT*log(γx[i])
                g_ideal_x = x[i]*RT*log(x[i])
                g_pure_x = x[i]*g_pures[i]
                gibbs += (g_E_x + g_ideal_x + g_pure_x)*(1-β)/RT
            end
        end
    end

    if is_vapour(phasey)
        lnϕy, voly = lnϕ(gas_model(fluidmodel), p, T, y, cache; phase=phasey, vol0=voly)
        for i in eachindex(lnK)
            if non_inx[i]
                lnK[i] = Inf
            elseif non_iny[i]
                lnK[i] = -Inf
            else
                ϕli = fug[i]
                p_i = sats[i][1]
                lnKi = log(γx[i]*p_i*ϕli/p) - lnϕy[i]
                !is_ideal && (lnKi += volx*(p - p_i)/RT) #add poynting corrections only if the fluid model itself has non-ideal corrections
                lnK[i] = lnKi
            end
            !non_iny[i] && β !== nothing && (gibbs += β*y[i]*(log(y[i]) + lnϕy[i]))
        end
    else
        γy = activity_coefficient(model, p, T, y)
        lnK .= log.(γx ./ γy)
        voly = volume(fluidmodel, p, T, y, phase = phasey, vol0 = voly)
        if β !== nothing
            for i in eachindex(y)
                if !non_inx[i]
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

function __tpflash_gibbs_reduced(wrapper::PTFlashWrapper{<:GammaPhi},p,T,x,y,β,eq,vols)
    pures = wrapper.model.fluid.pure
    model = wrapper.model
    fluidmodel = model.fluid.model
    g_pures = wrapper.μ
    RT = R̄*T

    gibbs = zero(Base.promote_eltype(model,p,T,x,β))

    if !isone(β)
        γx = activity_coefficient(model.activity, p, T, x)
        n = length(model)
        g_E_x = sum(x[i]*RT*log(γx[i]) for i ∈ 1:n)
        g_ideal_x = sum(x[i]*RT*(log(x[i])) for i ∈ 1:n)
        g_pure_x = dot(x,g_pures)
        gibbs += (g_E_x + g_ideal_x + g_pure_x)*(1-β)/RT
    end

    if is_vle(eq) && !iszero(β)
        gibbs += gibbs_free_energy(gas_model(fluidmodel),p,T,y,phase =:v)*β/R̄/T
    elseif !iszero(β) #lle
        γy = activity_coefficient(model.activity, p, T, y)
        g_E_y = sum(y[i]*RT*log(γy[i]) for i ∈ 1:n)
        g_ideal_y = sum(y[i]*R̄*T*(log(y[i])) for i ∈ 1:n)
        g_pure_y = dot(y,g_pures)
        gibbs += (g_E_y + g_ideal_y + g_pure_y)*β/RT
    end
    return gibbs
    #(gibbs_free_energy(model,p,T,x)*(1-β)+gibbs_free_energy(model,p,T,y)*β)/R̄/T
end

#TODO: derive expressions for this

function dgibbs_obj!(model::PTFlashWrapper{<:GammaPhi}, p, T, z, phasex, phasey,
    nx, ny, vcache, ny_var = nothing, in_equilibria = FillArrays.Fill(true,length(z)), non_inx = in_equilibria, non_iny = in_equilibria;
    F=nothing, G=nothing, H=nothing)
    throw(error("γ-ϕ Composite Model don't support gibbs energy optimization in MichelsenTPFlash."))
end

function K0_lle_init(wrapper::PTFlashWrapper{<:GammaPhi},p,T,z)
    return K0_lle_init(wrapper.model.activity,p,T,z)
end

function __eval_G_DETPFlash(wrapper::PTFlashWrapper{<:GammaPhi},p,T,x,equilibrium)
    model = wrapper.model
    phase = is_lle(equilibrium) ? :liquid : :unknown
    n = length(model)
    g_pures = wrapper.μ
    R = Rgas()
    RT = R*T
    γx = activity_coefficient(model.activity, p, T, x)
    g_E_x = sum(x[i]*RT*log(γx[i]) for i ∈ 1:n)
    g_ideal_x = sum(x[i]*RT*(log(x[i])) for i ∈ 1:n)
    g_pure_x = dot(x,g_pures)
    gl = (g_E_x + g_ideal_x + g_pure_x)
    vl = volume(model.fluid.model,p,T,x,phase = :l)
    if phase == :liquid
        return gl,vl
    else
        throw(error("γ-ϕ Composite Model does not support VLE calculation with `DETPFlash`. if you want to calculate LLE equilibria, try using `DETPFlash(equilibrium = :lle)`"))
        #=
        vv = volume(model.fluid.model,p,T,x,phase = :v)
        gv = VT_gibbs_free_energy(model.fluid.model, vv, T, x)
        if gv > gl
            return gl,vl
        else
            return gv,vv
        end
        =#
    end
end

#=
TPD support.

TODO: support vle in TPD.
=#
function tpd_delta_d_vapour!(d,wrapper,p,T)
    gas_model(wrapper) isa IdealModel && return d
    ϕsat,sat = wrapper.fug,wrapper.sat
    RT = R̄*T
    for i in eachindex(d)
        ps,vl,vv = sat[i]
        d[i] = d[i] - vl*(p - ps)/RT - log(ϕsat[i]) - log(ps/p)
    end
    return d
end

function tpd_input_composition(wrapper::PTFlashWrapper{<:GammaPhi},p,T,z,lle,cache = tpd_cache(model,p,T,z,di))

    TT = Base.promote_eltype(wrapper.model,p,T,z)

    pures = wrapper.model.fluid.pure
    model = wrapper.model
    fluidmodel = model.fluid.model
    g_pures = wrapper.μ
    RT = R̄*T

    d_l,d_v,_,_,_,Hϕ = cache

    TT = Base.promote_eltype(model,p,T,z)

    n = sum(z)
    logsumz = log(n)
    d,vl = _tpd_fz_and_v!(last(cache),wrapper,p,T,z,nothing,false,:liquid)
    d_l .= d
    d_l .+= log.(z) .- logsumz

    lle && return copy(d_l),:liquid,vl

    d,vv = _tpd_fz_and_v!(last(cache),wrapper,p,T,z,nothing,false,:vapour)
    d_v .= d
    d_v .+= log.(z) .- logsumz
    gr_l = dot(z,d_l)
    gr_v = dot(z,d_v)
    if gr_l < gr_v
        return copy(d_l),:liquid,vl
    else
        return copy(d_v),:vapour,vv
    end
end

function _tpd_fz_and_v!(cache,wrapper::PTFlashWrapper{<:GammaPhi},p,T,w,vol0,liquid_overpressure = false,phase = :l,_vol = nothing)
    pures = wrapper.model.fluid.pure
    model = wrapper.model
    fluidmodel = model.fluid.model
    g_pures = wrapper.μ
    RT = R̄*T

    if is_liquid(phase)
        γ = activity_coefficient(model.activity,p,T,w)
        v = one(eltype(γ))
        γ .= log.(γ)
        return γ,v,true
    else
        if _vol != nothing
            vol = _vol
        else
            vol = volume(gas_model(wrapper),p,T,w,phase = :vapour)
        end
        fxy,v = lnϕ!(cache,gas_model(wrapper),p,T,w;vol)
        overpressure = false
        if isnan(v) && liquid_overpressure && is_liquid(phase)
            overpressure = true
            #michelsen recomendation: when doing tpd, sometimes, the liquid cannot be created at the
            #specified conditions. try elevating the pressure at the first iter.
            fxy,v = lnϕ!(cache,gas_model(wrapper),1.2p,T,w,phase = phase)
        end
        tpd_delta_d_vapour!(fxy,wrapper,p,T)
        return fxy,v,overpressure
    end
end

function tpd_obj(wrapper::PTFlashWrapper{<:GammaPhi}, p, T, di, phasew, cache)

    function f(α)
        w,dtpd,_,_,vcache,Hϕ = cache
        ϕsat,sat = wrapper.fug,wrapper.sat
        w .= α .* α .* 0.25
        w ./= sum(w)
        if is_liquid(phasew)
            logγ,_ = _tpd_fz_and_v!(last(cache),wrapper,p,T,w,nothing,false,:liquid)
            dtpd .= log.(w) .+ logγ .- di
        else
            volw0 = vcache[]
            lnϕw, volw = lnϕ!(Hϕ, gas_model(wrapper), p, T, w; phase=phasew, vol0=volw0)
            tpd_delta_d_vapour!(lnϕw,wrapper,p,T)
            dtpd .= log.(w) .+ lnϕw .- di
            vcache[] = volw
        end
        fx = dot(w,dtpd) - sum(w) + 1
    end

    function g(df,α)
        w,dtpd,_,_,vcache,Hϕ = cache
        w .= α .* α .* 0.25
        w ./= sum(w)
        #@show w
        if is_liquid(phasew)
            logγ,_ = _tpd_fz_and_v!(last(cache),wrapper,p,T,w,nothing,false,:liquid)
            dtpd .= log.(w) .+ logγ .- di
        else
            volw0 = vcache[]
            lnϕw, volw = lnϕ!(Hϕ, gas_model(wrapper), p, T, w; phase=phasew, vol0=volw0)
            tpd_delta_d_vapour!(lnϕw,wrapper,p,T)
            dtpd .= log.(w) .+ lnϕw .- di
            vcache[] = volw
        end
        df .= dtpd .*  sqrt.(w)
        return df
    end

    function fg(df,α)
        w,dtpd,_,_,vcache,Hϕ = cache
        w .= α .* α .* 0.25
        w ./= sum(w)
        if is_liquid(phasew)
            logγ,_ = _tpd_fz_and_v!(last(cache),wrapper,p,T,w,nothing,false,:liquid)
            dtpd .= log.(w) .+ logγ .- di
        else
            volw0 = vcache[]
            lnϕw, volw = lnϕ!(Hϕ, gas_model(wrapper), p, T, w; phase=phasew, vol0=volw0)
            tpd_delta_d_vapour!(lnϕw,wrapper,p,T)
            dtpd .= log.(w) .+ lnϕw .- di
            vcache[] = volw
        end
        df .= dtpd .*  sqrt.(w)
        fx = dot(w,dtpd) - sum(w) + 1
        return fx,df
    end

    #=
    from thermopack:
    We see that ln Wi + lnφ(W) − di will be zero at the solution of the tangent plane minimisation.
    It can therefore be removed from the second derivative, without affecting the convergence properties.

    function fgh(df,d2f,α)
        w,dtpd,_,_,vcache,Hϕ = cache
        nc = length(model)
        w .= α .* α .* 0.25
        w ./= sum(w)
        volw0 = vcache[]
        lnϕw, ∂lnϕ∂nw, ∂lnϕ∂Pw, volw = ∂lnϕ∂n∂P(model, p, T, w, Hϕ; phase=phase, vol0=volw0)
        for i in 1:nc
            xi = sqrt(w[i])
            for j in 1:nc
                xj = sqrt(w[j])
                δij = Int(i == j)
                d2f[i,j] = δij + xi*xj*∂lnϕ∂nw[i,j]
            end
        end
        dtpd .= log.(w) .+ lnϕw .- di
        df .= dtpd .*  sqrt.(w)
        fx = dot(w,dtpd) - sum(w) + 1
        #if fx < -1e-10 && break_first
        #    df .= 0
        #end
        vcache[] = volw
        return fx,df,d2f
    end

    function h(d2f,α)
        w,dtpd,_,_,vcache,Hϕ = cache
        nc = length(model)
        w .= α .* α .* 0.25
        w ./= sum(w)
        volw0 = vcache[]
        lnϕw, ∂lnϕ∂nw, ∂lnϕ∂Pw, volw = ∂lnϕ∂n∂P(model, p, T, w, Hϕ; phase=phase, vol0=volw0)
        for i in 1:nc
            xi = sqrt(w[i])
            for j in 1:nc
                xj = sqrt(w[j])
                δij = Int(i == j)
                d2f[i,j] = δij + xi*xj*∂lnϕ∂nw[i,j]# + 0.5*αi*dtpd[i]
            end
        end
        vcache[] = volw
        d2f
    end
    =#
    obj = NLSolvers.ScalarObjective(f=f,g=g,fg=fg,fgh=nothing,h=nothing)
    optprob = OptimizationProblem(obj = obj,inplace = true)
end

function tpd_optimization(model::PTFlashWrapper{<:GammaPhi},p,T,z,w0,di,cache = tpd_cache(model,p,T,z,K0),phasew = :liquid)
    w,_,_,dzz,vcache,Hϕ = cache
    α0 = 2 .* sqrt.(w0)
    prob = tpd_obj(model, p, T, dzz, phasew, cache)
    opt_options = OptimizationOptions(f_abstol = 1e-12,f_reltol = 1e-10,maxiter = 100)
    lb,ub = similar(α0),similar(α0)
    lb .= 0
    ub .= Inf
    res = Solvers.optimize(prob, α0, Solvers.LineSearch(Solvers.BFGS(),Solvers.BoundedLineSearch(lb,ub)), opt_options)
    α = Solvers.x_sol(res)
    w .= α .* α .* 0.25
    w ./= sum(w)
    tpd = Solvers.x_minimum(res)
    vw = vcache[]
    return w,tpd,vw
end

export GammaPhi
