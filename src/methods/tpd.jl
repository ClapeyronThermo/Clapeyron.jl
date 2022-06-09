
function tpd_obj!(model::EoSModel, p, T, di, α, phasew; volw0=nothing,
                  F=nothing, G=nothing, H=nothing)
    # Function that computes the TPD function, its gradient and its hessian
    nc = length(model)
    w = α.^2 /4.
    if H !== nothing
        # computing the hessian
        lnϕw, ∂lnϕ∂nw, ∂lnϕ∂Pw, volw = ∂lnϕ∂n∂P(model, p, T, w; phase=phasew, vol0=volw0)
        if isnan(volw)
            lnϕw, ∂lnϕ∂nw, ∂lnϕ∂Pw, volw = ∂lnϕ∂n∂P(model, p, T, w; phase=phasew, vol0=nothing)
        end
        dtpd = log.(w) + lnϕw - di
        gi = dtpd.*(α./2)
        eye = LinearAlgebra.I(nc)
        #TODO: check that just using Identity without instantiation works
        H .= eye .* (1. .+  (gi./α))  .+ sqrt.(w * w') .* ∂lnϕ∂nw
    else
        lnϕw, volw = lnϕ(model, p, T, w; phase=phasew, vol0=volw0)
        if isnan(volw)
            lnϕw, volw = lnϕ(model, p, T, w; phase=phasew, vol0=nothing)
        end
        dtpd = log.(w) + lnϕw - di
        gi = dtpd.*(α./2)
    end

    if G !== nothing
        # computing the gradient
         G .= gi
    end
    if F !== nothing
        # computing the TPD value
        #tpdi = w_notzero.*(dtpd .- 1.)
        #tpd = sum(tpdi) + 1
        tpd = dot(w,dtpd) - sum(w) + 1
        return tpd
    end
end

function tpd_ss(model::EoSModel, p, T, di, w, phasew; volw0=nothing,max_iters = 5)
    # Function that minimizes the tpd function by Successive Substitution
    volw0 === nothing && (volw0 = volume(model, p, T, w, phase = phasew))
    volw = volw0
    wres = copy(w)
    lnw = copy(w)
    for _ in 1:max_iters
        lnϕw, volw = lnϕ(model, p, T, wres; phase=phasew, vol0=volw)
        lnw .= di .- lnϕw
        wres .= exp.(lnw)
        wres ./= sum(w)
    end
    return w, volw
end

function tpd_min(model::EoSModel, p, T, z, w, phasez, phasew; volz0=nothing, volw0=nothing)
    # Function that minimizes the tpd function first by Successive Substitution
    # and then by a Newton's method
    # out = minimized trial phase composition (w) and its tpd value
    nc = length(model)
    # computing the di vector for the phase z (constant along the minimization)
    lnϕz, volz = lnϕ(model, p, T, z; phase=phasez, vol0=volz0)
    di = log.(z) + lnϕz

    volw = volw0
    # improving initial guess by Successive Substitution
    w, volw = tpd_ss(model, p, T, di, w, phasew; volw0=nothing)

    # change of variable to "number of moles"
    α0 = max.(2 .* sqrt.(w),one(eltype(w))*1e-2)

    # minimizing the TPD by Newton's method
    dftpd!(F, G, H, α) = tpd_obj!(model, p, T, di, α, phasew, volw0=volw, F=F, G=G, H=H)
    sol = Solvers.optimize(Solvers.only_fgh!(dftpd!), α0) #, method=LineSearch(Newton()))#  , Optim.Newton())
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

    #TODO for the future:
    #this operation is a "embarrasingly parallel" problem, multithreading will surely speed this up
    #but Base.@threads on julia 1.6 behaves on a static manner, on 1.8 onwards, there is Base.@threads :dynamic,
    #that allows nesting. Ideally, all Clapeyron operations should be multithread-friendly.
    for (phasez,phasew) in phasepairs
        for i in 1:length(model)
            w0 = Id[i, :]
            try
                w, tpd = tpd_min(model,p,T,z,w0,phasez,phasew)
                if tpd < 0. && ~isapprox(z, w, atol=1e-3)
                    if length(w_array) > 0
                        already_computed = false
                        for ws in w_array
                            # check if the minimum is already stored
                            already_computed = already_computed || isapprox(ws, w, atol=1e-3)
                        end
                        if ~already_computed
                            push!(w_array, index_expansion(w,z_notzero))
                            push!(tpd_array, tpd)
                            push!(phasez_array, phasez)
                            push!(phasew_array, phasew)
                            # println(i, ' ', phasez, ' ', phasew, ' ', w, ' ', tpd)
                        end
                    else
                        push!(w_array, w)
                        push!(tpd_array, tpd)
                        push!(phasez_array, phasez)
                        push!(phasew_array, phasew)
                        # println(i, ' ', phasez, ' ', phasew, ' ', w, ' ', tpd)
                    end
                end
            catch e
                # If the minimization is not successful print this line
                verbose && @warn("""Failed to minimize the TPD function for w0 = $(Id[i, :]),
                        phasew = $phasew, phasez = $phasez at pressure p = $p [Pa],
                        T = $T [K] and global composition = $z""")
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
    tpd(model,p,T,z;verbose=false)

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

export tpd