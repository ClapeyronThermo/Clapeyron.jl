


function tpd_obj!(model::EoSModel, p, T, di, α, phasew, z_notzero; volw0=nothing,
                  F=nothing, G=nothing, H=nothing)
    # Function that computes the TPD function, its gradient and its hessian
    nc = length(model)
    ncomponents = length(α)

    w = zeros(nc)
    w[z_notzero] = α.^2 /4.
    w_notzero = w[z_notzero]

    if H !== nothing
        # computing the hessian
        lnϕw, ∂lnϕ∂nw, ∂lnϕ∂Pw, volw = ∂lnϕ∂n∂P(model, p, T, w; phase=phasew, vol0=volw0)
        if isnan(volw)
            lnϕw, ∂lnϕ∂nw, ∂lnϕ∂Pw, volw = ∂lnϕ∂n∂P(model, p, T, w; phase=phasew, vol0=nothing)
        end

        dtpd = log.(w_notzero) + lnϕw[z_notzero] - di
        gi = dtpd.*(α./2)

        eye = Identity(ncomponents)
        #TODO: check that just using Identity without instantiation works
        hess = @view H[:, :]
        hess .= eye .* (1. .+  (gi./α))  .+ sqrt.(w_notzero * w_notzero') .* ∂lnϕ∂nw[z_notzero, z_notzero]
    else
        lnϕw, volw = lnϕ(model, p, T, w; phase=phasew, vol0=volw0)
        if isnan(volw)
            lnϕw, volw = lnϕ(model, p, T, w; phase=phasew, vol0=nothing)
        end
        dtpd = log.(w_notzero) + lnϕw[z_notzero] - di
        gi = dtpd.*(α./2)
    end

    if G !== nothing
        Gvec = vec(G)
        # computing the gradient
        Gvec .= gi
    end

    if F !== nothing
        # computing the TPD value
        #tpdi = w_notzero.*(dtpd .- 1.)
        #tpd = sum(tpdi) + 1
        tpd = dot(w_notzero,dtpd) - sum(w_notzero) + 1
        return tpd
    end
end

function tpd_ss(model::EoSModel, p, T, di, w, phasew, z_notzero; volw0=nothing)
    # Function that minimizes the tpd function by Successive Substitution
    nc = length(model)
    volw = volw0
    for i in 1:5
        lnϕw, volw = lnϕ(model, p, T, w; phase=phasew, vol0=volw)
        lnw = di - lnϕw[z_notzero]
        w = zeros(nc)
        w_notzero = exp.(lnw)
        w[z_notzero] = exp.(lnw)
        w /= sum(w)
    end
    return w, volw
end

function tpd_min(model::EoSModel, p, T, z, w, phasez, phasew; volz0=nothing, volw0=nothing)
    # Function that minimizes the tpd function first by Successive Substitution
    # and then by a Newton's method
    # out = minimized trial phase composition (w) and its tpd value

    z_notzero = z .> 0.
    # computing the di vector for the phase z (constant along the minimization)
    lnϕz, volz = lnϕ(model, p, T, z; phase=phasez, vol0=volz0)
    di = log.(z[z_notzero]) + lnϕz[z_notzero]

    volw = volw0
    # improving initial guess by Successive Substitution
    w, volw = tpd_ss(model, p, T, di, w, phasew, z_notzero; volw0=nothing)

    # change of variable to "number of moles"
    α0 = 2*w.^0.5
    α0[α0 .< 1e-2] .= 1e-2
    α0 = α0[z_notzero]

    # minimizing the TPD by Newton's method
    dftpd!(F, G, H, α) = tpd_obj!(model, p, T, di, α, phasew, z_notzero, volw0=volw, F=F, G=G, H=H)
    # sol = Optim.optimize(only_fgh!(dftpd!), α0, Optim.Newton())
    sol = Solvers.optimize(Solvers.only_fgh!(dftpd!), α0) #, method=LineSearch(Newton()))#  , Optim.Newton())
    nc = length(model)
    # computing phase composition
    w = zeros(nc)
    w[z_notzero] = sol.info.solution.^2 / 4.
    w = w/sum(w)
    tpd = sol.info.minimum
    return w, tpd
end

function all_tpd(model::EoSModel, p, T, z)
    # Function that minimizes the tpd function first by Successive Substitution
    # and then by a Newton's method multiple times
    # This functions attempts to find all the tpd minima
    # out = array of minimas composition (w), array of tpd values (tpd), array of phase z state, array of phase w state

    z_notzero = z .> 0.

    nc = length(model)
    Id = Matrix{Float64}(Identity, nc, nc)

    w_array = []
    tpd_array = []
    phasez_array = []
    phasew_array = []

    for phasez in (:liquid, :vapor)
        # computing the di vector for the phase z (constant along the minimization for a given phasez)
        lnϕz, volz = lnϕ(model, p, T, z; phase=phasez, vol0=nothing)
        di = log.(z[z_notzero]) + lnϕz[z_notzero]

        if phasez == :liquid
            possible_phasew = [:liquid, :vapor]
        else
            possible_phasew = [:liquid]
        end

        for phasew in possible_phasew
            for i in 1:nc
                volw = nothing
                w = Id[i, :]
                # improving initial guess by Successive Substitution
                w, volw = tpd_ss(model, p, T, di, w, phasew, z_notzero; volw0=nothing)

                # change of variable to "number of moles"
                α0 = 2*w.^0.5
                α0[α0 .< 1e-2] .= 1e-2
                α0 = α0[z_notzero]
                try
                    # minimizing the TPD by Newton's method
                    dftpd!(F, G, H, α) = tpd_obj!(model, p, T, di, α, phasew, z_notzero, volw0=volw, F=F, G=G, H=H)
                    # sol = Optim.optimize(only_fgh!(dftpd!), α0, Optim.Newton())
                    sol = Solvers.optimize(Solvers.only_fgh!(dftpd!), α0)
                    w = zeros(nc)
                    w[z_notzero] = sol.info.solution.^2 / 4.
                    w = w/sum(w)
                    tpd = sol.info.minimum
                    if tpd < 0. && ~isapprox(z, w, atol=1e-3)
                        if length(w_array) > 0
                            already_computed = false
                            for ws in w_array
                                # check if the minimum is already stored
                                already_computed = already_computed || isapprox(ws, w, atol=1e-3)
                            end
                            if ~already_computed
                                push!(w_array, w)
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
                catch
                    # If the minimization is not successful print this line
                    println("""Failed to minimize the TPD function for w0 = $(Id[i, :]),
                            phasew = $phasew, phasez = $phasez at pressure p = $p [Pa],
                            T = $T [K] and global composition = $z""")
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

function lle_init(model::EoSModel, p, T, z)
    # Function that minimizes the tpd function first by Successive Substitution
    # and then by a Newton's method multiple times
    # This functions attempts to find all the liquid tpd minima
    # out = array of minimas composition, array of tpd values
    _1 = one(p+T+first(z))
    TYPE = typeof(_1)
    z_notzero = z .> 0.

    nc = length(model)
    Id = Matrix{TYPE}(Identity, nc, nc)

    w_array = []
    tpd_array = []

    phasez = :liquid
    phasew = :liquid
    # computing the di vector for the phase z (constant along the minimization for a given phasez)
    lnϕz, volz = lnϕ(model, p, T, z; phase=phasez, vol0=nothing)
    di = log.(z[z_notzero]) + lnϕz[z_notzero]

    for i in 1:nc
        volw = nothing
        w = Id[i, :]
        # improving initial guess by Successive Substitution
        w, volw = tpd_ss(model, p, T, di, w, phasew, z_notzero; volw0=nothing)

        α0 = 2*w.^0.5
        α0[α0 .< 1e-2] .= 1e-2
        α0 = α0[z_notzero]
        try
            # minimizing the TPD by Newton's method
            dftpd!(F, G, H, α) = tpd_obj!(model, p, T, di, α, phasew, z_notzero, volw0=volw, F=F, G=G, H=H)
            # sol = Optim.optimize(only_fgh!(dftpd!), α0, Optim.Newton())
            sol = Solvers.optimize(Solvers.only_fgh!(dftpd!), α0)
            w = zeros(nc)
            w[z_notzero] = sol.info.solution.^2 / 4.
            w = w/sum(w)
            tpd = sol.info.minimum
            if tpd < 0. && ~isapprox(z, w, atol=1e-3)
                if length(w_array) > 0
                    already_computed = false
                    for ws in w_array
                        # check if the minimum is already stored
                        already_computed = already_computed || isapprox(ws, w, atol=1e-3)
                    end
                    if ~already_computed
                        push!(w_array, w)
                        push!(tpd_array, tpd)
                    end

                else
                    push!(w_array, w)
                    push!(tpd_array, tpd)
                end

            end
        catch
            # If the minimization is not successful print this line
            println("""Failed to minimize the TPD function for w0 = $(Id[i, :]),
                    phasew = $phasew, phasez = $phasez at pressure p = $p [Pa],
                    T = $T [K] and global composition = $z""")
        end
    end

    index = sortperm(tpd_array)
    w_array = w_array[index]
    tpd_array = tpd_array[index]

    return w_array, tpd_array
end
