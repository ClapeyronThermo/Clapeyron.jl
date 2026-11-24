const NO_FLUID_SOLID_MSG = "MultiPhaseTPFlash: no fluid phase in equilibrium with the solids, use CEMFlash"

function tp_flash_impl(model::CompositeModel{<:Any,<:EoSModel},p,T,z,method::MultiPhaseTPFlash)
    return tp_flash_multi_solid(model,p,T,z,method,Int[],Float64[])
end

#=
Pure (stoichiometric) solid phases. The amounts s of the active solids (formula units per mole of feed) are solved with a
damped Newton method on the solid equilibrium conditions Σₖνₖ(ln fₖ - ln fₖ⁰) = Δμ/RT (pure liquid reference, sle_Δμ),
with an ideal solution jacobian. Each evaluation solves the fluid phases with method on the remaining feed z - νs.
Solids are added (most supersaturated first) and removed (amount reaches zero) until all solids are stable.
=#
function tp_flash_multi_solid(model::CompositeModel,p,T,z,method,active,s)
    fluid = __tpflash_cache_model(model.fluid,p,T,z,:unknown)
    N,ns = length(model.fluid),length(model.solid)
    ν = zeros(N,ns)
    for j in 1:ns
        idx,νj = sle_mapping(model,j)
        ν[idx,j] .= νj
    end
    Δμ = [sle_Δμ(model,p,T,j) for j in 1:ns]
    cache = ∂lnϕ_cache(fluid,p,T,z,Val{false}())
    insolid = [any(>(0),@view(ν[k,:])) for k in 1:N]
    lnf0 = [insolid[k] ? first(modified_lnϕ(fluid,p,T,Float64.(1:N .== k),cache;phase = :liquid))[k] : 0.0 for k in 1:N]
    function fluid_state(s)
        nf = z .- @view(ν[:,active])*s
        res = tp_flash_impl(model.fluid,p,T,nf ./ sum(nf),method)
        x = res.compositions[1]
        lnϕ,_ = modified_lnϕ(fluid,p,T,x,cache;phase = identify_phase(res,1),vol0 = res.volumes[1])
        return nf,res,lnϕ .+ log.(x)
    end
    #Δμ/RT - Σₖνₖ(ln fₖ - ln fₖ⁰): negative if the solid is supersaturated
    D(lnf,j) = Δμ[j] - sum(ν[k,j]*(lnf[k] - lnf0[k]) for k in 1:N if ν[k,j] > 0)
    tol = method.K_tol
    nf,res,lnf = fluid_state(s)
    converged = false
    for _ in 1:method.phase_iters
        for _ in 1:100
            F = Float64[-D(lnf,j) for j in active]
            norm(F,Inf) < tol && break
            Nf = sum(nf)
            J = Float64[sum(ν[k,i]*(sum(@view(ν[:,j]))/Nf - ν[k,j]/nf[k]) for k in 1:N if ν[k,i] > 0) for i in active, j in active]
            Jf = lu!(J,check = false)
            issuccess(Jf) || throw(ArgumentError(NO_FLUID_SOLID_MSG))
            Δs = -(Jf\F)
            νΔs = @view(ν[:,active])*Δs
            #a solid whose amount reaches zero is removed, the fluid feed stays positive
            αs,iremove = 1.0,0
            for i in eachindex(Δs)
                Δs[i] < 0 && -s[i]/Δs[i] < αs && ((αs,iremove) = (-s[i]/Δs[i],i))
            end
            α = Solvers.positive_linesearch(nf,-νΔs,τ = 0.99,decay = 1.0)
            if iremove != 0 && αs <= α
                s .+= αs .* Δs
                deleteat!(active,iremove)
                deleteat!(s,iremove)
                nf,res,lnf = fluid_state(s)
                continue
            end
            normF = norm(F)
            for _ in 1:20
                s_new = s .+ α .* Δs
                nf_new,res_new,lnf_new = fluid_state(s_new)
                if norm(Float64[D(lnf_new,j) for j in active]) < normF || α < 1e-4
                    s .= s_new
                    nf,res,lnf = nf_new,res_new,lnf_new
                    break
                end
                α /= 2
            end
        end
        if !all(j -> abs(D(lnf,j)) < tol,active)
            sum(nf) < 1e-6 && throw(ArgumentError(NO_FLUID_SOLID_MSG))
            break
        end
        inactive = setdiff(1:ns,active)
        isempty(inactive) && (converged = true; break)
        Dmin,i = findmin(j -> D(lnf,j),inactive)
        Dmin > -tol && (converged = true; break)
        #phase rule: a fluid phase coexists with at most N - 1 solids (except at invariant points)
        length(active) + 1 >= N && throw(ArgumentError(NO_FLUID_SOLID_MSG))
        push!(active,inactive[i])
        push!(s,0.0)
    end
    converged || error("MultiPhaseTPFlash: the solid-fluid equilibrium did not converge")
    Nf = sum(nf)
    #the fluid flash can return copies of a stable phase
    idx = Int[]
    β = eltype(res.fractions)[]
    for i in 1:numphases(res)
        l = findfirst(m -> identify_phase(res,m) == identify_phase(res,i) && isapprox(res.compositions[m],res.compositions[i],rtol = 1e-8),idx)
        l === nothing ? (push!(idx,i);push!(β,res.fractions[i])) : (β[l] += res.fractions[i])
    end
    comps = [collect(res.compositions[i]) for i in idx]
    β .*= Nf/sum(β)
    vols = res.volumes[idx]
    vapour_idx = res.data.vapour_idx > 0 ? findfirst(==(res.data.vapour_idx),idx) : res.data.vapour_idx
    g = res.data.g*Nf
    gref = [insolid[k] ? first(modified_gibbs(fluid,p,T,Float64.(1:N .== k),:liquid)) : 0.0 for k in 1:N]
    for (j,sj) in zip(active,s)
        νj = ν[:,j]
        ∑ν = sum(νj)
        push!(comps,νj ./ ∑ν)
        push!(β,sj*∑ν)
        push!(vols,volume(sle_pure_solid(model.solid,j),p,T,SA[1.0])/∑ν)
        g += sj*(dot(νj,gref) + Δμ[j])
    end
    solid_idx = isempty(active) ? -1 : length(idx) + 1
    return FlashResult(comps,β,vols,FlashData(p,T,g,vapour_idx,solid_idx))
end
