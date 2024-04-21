
function tpd_obj!(model::EoSModel, p, T, di, α, phasew, vcache;
                  F=nothing, G=nothing, H=nothing)
    # Function that computes the TPD function, its gradient and its hessian
    nc = length(model)
    w = α.^2 /4.0
    #sqrt(w) = 0.5*α
    volw0 = vcache[]
    if H !== nothing
        # computing the hessian
        lnϕw, ∂lnϕ∂nw, ∂lnϕ∂Pw, volw = ∂lnϕ∂n∂P(model, p, T, w; phase=phasew, vol0=volw0)
        for i in 1:nc
            αi = α[i]
            for j in 1:nc
                αj = α[j]
                δij = Int(i == j)

                #=
                from thermopack:
                We see that ln Wi + lnφ(W) − di will be zero at the solution of the tangent plane minimisation.
                It can therefore be removed from the second derivative, without affecting the convergence properties.
                =#

                H[i,j] += δij + 0.25*αi*αj*∂lnϕ∂nw[i,j] #+ 0.5*αi*dtpd[i]
            end
        end
    else
        lnϕw, volw = lnϕ(model, p, T, w; phase=phasew, vol0=volw0)
    end
    dtpd = log.(w) + lnϕw - di
    vcache[] = volw
    if G !== nothing
        # computing the gradient
         G .= dtpd.* ( α./ 2)
    end
    if F !== nothing
        # computing the TPD value
        #tpdi = w_notzero.*(dtpd .- 1.)
        #tpd = sum(tpdi) + 1
        tpd = dot(w,dtpd) - sum(w) + 1
        return tpd
    end
end

function tpd_ss(model::EoSModel, p, T, di, w0, phasew; volw0=nothing,max_iters = 10)
    # Function that minimizes the tpd function by Successive Substitution
    volw0 === nothing && (volw0 = volume(model, p, T, w0, phase = phasew))
    volw = volw0
    w = copy(w0)
    lnw = copy(w0)
    tpd = zero(T+p+first(w0))
    for _ in 1:max_iters
        lnϕw, volw = lnϕ(model, p, T, w; phase=phasew, vol0=volw)
        lnw .= di .- lnϕw
        w .= exp.(lnw)
        w ./= sum(w)
        !isfinite(first(w)) && break
    end
    return w, volw
end

function tpd_min(model::EoSModel, p, T, di, z, w0, phasez, phasew; volz0=nothing, volw0=nothing)
    # Function that minimizes the tpd function first by Successive Substitution
    # and then by a Newton's method
    # out = minimized trial phase composition (w) and its tpd value
    nc = length(model)
    volw = volw0
    # improving initial guess by Successive Substitution
    w, volw = tpd_ss(model, p, T, di, w0, phasew; volw0=nothing)

    if !isfinite(volw) || !isfinite(first(w)) #iteration returned non finite values
        return w, zero(eltype(w))/zero(eltype(w))
    end

    if isapprox(z, w, atol=1e-3) #ss iteration converged to z
        return w, zero(eltype(w))/zero(eltype(w))
    end

    # change of variable to "number of moles"
    #α0 = max.(2 .* sqrt.(w),one(eltype(w))*1e-8)
    α0 = 2 .* sqrt.(w)
    vcache = Ref(volw)
    # minimizing the TPD by Newton's method
    dftpd!(F, G, H, α) = tpd_obj!(model, p, T, di, α, phasew, vcache, F=F, G=G, H=H)
    sol = Solvers.optimize(Solvers.only_fgh!(dftpd!), α0,LineSearch(Newton())) #, method=LineSearch(Newton()))#  , Optim.Newton())
    # computing phase composition
    w = sol.info.solution.^2 / 4
    w ./= sum(w)
    tpd = sol.info.minimum
    return w, tpd
end

function all_tpd(model::EoSModel, p, T, z,phasepairs = ((:liquid,:vapour),(:liquid,:liquid),(:vapour,:liquid));verbose = false)
    # Function that minimizes the tpd function first by Successive Substitution
    model_full,z_full = model,z
    model, z_notzero = index_reduction(model_full,z_full)
    z = z_full[z_notzero]
    nc = length(model_full)

    _1 = one(p+T+first(z))
    nc = length(model)
    Id = fill(zero(eltype(z)),nc,nc)
    for i in diagind(Id)
        Id[i] = 1.0
    end
    w_array = Vector{Vector{eltype(_1)}}(undef,0)
    tpd_array = fill(_1,0)
    phasez_array = fill(:x,0)
    phasew_array = fill(:x,0)

    #cache di
    di_dict = Dict{Symbol,Vector{Float64}}()
    #g_dict = Dict{Symbol,Float64}()
    for (phasez,_) in phasepairs
        if !haskey(di_dict,phasez)
            lnϕz, volz = lnϕ(model, p, T, z; phase=phasez)
            #g_dict[phasez] = VT_gibbs_free_energy(model,volz,T,z)
            isnan(volz) && continue
            di = log.(z) + lnϕz
            add_to = true
            for (phaseij,dij) in pairs(di_dict)
                if isapprox(di,dij,rtol = 1e-5) #new similar phases
                    add_to = false
                end
            end
            if add_to
                di_dict[phasez] = di
            end
        end
    end

    #TODO for the future:
    #this operation is a "embarrasingly parallel" problem, multithreading will surely speed this up
    #but Base.@threads on julia 1.6 behaves on a static manner, on 1.8 onwards, there is Base.@threads :dynamic,
    #that allows nesting. Ideally, all Clapeyron operations should be multithread-friendly.
    for (phasez,phasew) in phasepairs
        !haskey(di_dict,phasez) && continue
        di = di_dict[phasez]
        for i in 1:length(model)
            w0 = Id[i, :] #vector of single component
            w, tpd = tpd_min(model,p,T,di,z,w0,phasez,phasew)
            isnan(tpd) && continue
            if tpd < 0. && !isapprox(z, w, atol=1e-3)
                if length(w_array) == 0
                    push!(w_array, w)
                    push!(tpd_array, tpd)
                    push!(phasez_array, phasez)
                    push!(phasew_array, phasew)
                    continue
                end
                already_computed = false
                for ws in w_array
                    # check if the minimum is already stored
                    already_computed = already_computed || isapprox(ws, w, atol=1e-3)
                end
                if !already_computed
                    push!(w_array, index_expansion(w,z_notzero))
                    push!(tpd_array, tpd)
                    push!(phasez_array, phasez)
                    push!(phasew_array, phasew)
                    # println(i, ' ', phasez, ' ', phasew, ' ', w, ' ', tpd)
                end
            end
        end
    end
    # sort the obtained tpd minimas
    index = sortperm(tpd_array)
    w_array = w_array[index]
    tpd_array = tpd_array[index]
    phasez_array = phasez_array[index]
    phasew_array = phasew_array[index]
    return w_array, tpd_array, phasez_array, phasew_array
end

"""
    tpd(model,p,T,z;verbose = false)

Calculates the Tangent plane distance function (`tpd`). It returns:

- a vector with trial phase compositions where `tpd < 0`
- a vector with the `tpd` values
- a vector with symbols indicating the phase of the input composition
- a vector with symbols indicating the phase of the trial composition

It iterates over each two-phase combination, starting from pure trial compositions, it does succesive substitution, then Gibbs optimization.

If the vectors are empty, then the procedure couldn't find a negative `tpd`. That is an indication that the phase is (almost) surely stable.

"""
tpd(model,p,T,z;verbose = false) = all_tpd(model,p,T,z;verbose = verbose)

function lle_init(model::EoSModel, p, T, z;verbose = false)
    w_array, tpd_array, _, _ = all_tpd(model,p,T,z,((:liquid,:liquid),);verbose = verbose)
    return w_array, tpd_array
end

# function K0_lle_init(model::EoSModel, p, T, z)
#     w,_ = lle_init(model, p, T, z)
#     #vlle, or other things
#     if length(w) != 2
#         _0 = zero(eltype(w))
#         return fill(_0/_0,length(w))
#     else
#         x1,x2 = w[1],w[2]
#         return x1 ./ x2
#     end
# end

function K0_lle_init(model::EoSModel, p, T, z)
    nc = length(model)
    if nc == 2
        z_test = [1. 1e-6; 1e-6 1.]
    else
        z_test = ones(Int64(nc*(1+(nc-1)/2)),nc).*1e-3
        for i in 1:nc
            z_test[i,i] = 1.
        end
        k = nc+1
        for i in 1:nc-1
            for j in i+1:nc
                z_test[k,i] = 0.5
                z_test[k,j] = 0.5
                k += 1
            end
        end
    end

    z_test = z_test .* z'
    z_test = z_test ./ sum(z_test;dims=2)
    ntest = length(z_test[:,1])
    γ = zeros(ntest,nc)
    for i in 1:ntest
        γ[i,:] = activity_coefficient(model,p,T,z_test[i,:])
    end

    err = ones(ntest,ntest)*Inf
    for i in 1:ntest
        for j in i+1:ntest
            ϕi = (z_test[i,:].-z)./(z_test[i,:]-z_test[j,:])
            ϕ = sum(ϕi[isfinite.(ϕi)])/sum(isfinite.(ϕi))
            err[i,j] = sum(log.(γ[i,:].*z_test[i,:]) .-log.(γ[j,:].*z_test[j,:])+log.(γ[i,:].*z_test[i,:]).*(ϕ*z_test[i,:]+(1-ϕ)*z_test[j,:].-z))
        end
    end

    (val, idx) = findmin(err)

    K0 = γ[idx[1],:]./γ[idx[2],:]

    #=
    #extra step, reduce magnitudes if we aren't in a single phase
    g0 = dot(z, K0) - 1.
    g1 = 1. - sum(zi/Ki for (zi,Ki) in zip(z,K0))

    if g0 < 0 && g1 <= 0 #we need to correct the value of g0
        g0i = z .* K0
        kz,idx = findmax(g0i)
        #the maximum K value is too great
    elseif g1 > 0 && g0 >= 0 #we need to correct the value of g1
        g1i = z ./ K0

    else #bail out here?

    end =#
    return K0
end

#=
function K0_lle_init(model::EoSModel, p, T, z)
    w,_ = lle_init(model, p, T, z)
    #vlle, or other things
    if length(w) == 2
        x1,x2 = w[1],w[2]
        return x1 ./ x2
    elseif length(w) == 1
        x1 = w[1]
        x2 = Fractions.neg(x1)
        return x1 ./ x2
    else
        _0 = zero(eltype(w))
        return fill(_0/_0,length(w))
    end
end
=#

export tpd