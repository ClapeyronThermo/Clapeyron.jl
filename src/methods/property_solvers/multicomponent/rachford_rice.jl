"""
    rr_vle_vapor_fraction(K,z,α = NaN)

Given a vector of K values and a vector of compositions, calculates the vapor fraction `β`.
The algorithm is a modification of _(1)_ , with safeguards for extreme cases.

If the algorithm fails to converge, returns `NaN`. if it converges to a value `β ∉ [0,1]`, returns `-Inf` or `Inf`, depending on the case.

## References
1. Vassilis Gaganis, "Solution of the Rachford Rice equation using perturbation analysis", Fluid Phase Equilibria, Volume 536,2021,112981
"""
function rr_vle_vapor_fraction(K,z)
    #Solution of the Rachford Rice equation using perturbation analysis
    #https://doi.org/10.1016/j.fluid.2021.112981
    #strategy: find 4 strongest compositions
    #calculate exact β for those 4
    #real_β = β_ℵ1 + ε*β_ℵ2(ε), where ℵ1 is the strong group of component and
    #ℵ2 is the weak group of components.
    #on ε = 0, we have the aprox solution
    #on ε = 1 we have the real solution
    #ε*β_ℵ2(ε) ≈ ∑ci*ε^i
    #we obtain the coefficients and evaluate ε = 1.
    #for n = length(z) <= 4, we return exact solutions.
    n_z = length(z)
    _0 = zero(first(z)+first(K))
    _1 = one(_0)
    kmin = minimum(K)
    kmax = maximum(K)
    βmax = 1/(1-kmin)
    βmin = 1/(1-kmax)
    kmax < 1 && return (-_1/_0)
    kmin > 1 && return (_1/_0)
    if n_z <= 4 && all(Base.Fix2(>,0),z) #all z > 0
        return rr_vle_vapor_fraction_exact(K,z)
    end

    #z_ℵ1 is normalized
    ℵ1,k_ℵ1,z_ℵ1 = Clapeyron.rr_find_strongest(K,z)
    β_near0,β_near1 = extrema(((_1/(_1 -k_ℵ1[1])),(_1/(_1 -k_ℵ1[2]))))
    near_mean = sqrt(abs(β_near0)*abs(_1-β_near1))

    β0 = rr_vle_vapor_fraction_exact(k_ℵ1,z_ℵ1)
    b0 = β0
    M1,M2 = _0,_0
    N1,N2 = _0,_0
    Ξ1,Ξ2 = _0,_0
    P1 = _0
    Λ2 = _0
    for i in 1:length(z)
        βi = _1/(_1-K[i])
        γi = βi - b0
        γiinv = _1/γi
        λi = z[i]*γiinv
        μi = λi*γiinv
        νi = μi*γiinv
        ξi = νi*γiinv
        if i in ℵ1
            M1 += μi
            N1 += νi
            Ξ1 += ξi
            P1 += ξi*γiinv
        else
            M2 += μi
            N2 += νi
            Ξ2 += ξi
            Λ2 += λi
        end
    end
    b1 = -Λ2/M1
    b2 = -(b1*b1*N1 + b1*M2)/M1
    b11 = b1*b1
    b111 = b11*b1
    b12 = b1*b2
    b22 = b2*b2
    b3 = -((b111*Ξ1+2*b12*N1) + (b11*N2+b2*M2))/M1
    t1 = b11*b11*P1 + 3*b11*b2*Ξ1 + (b22+2*b1*b3)*N1
    t2 = b111*Ξ2 + 2*b12*N2 + b3*M2
    b4 = -(t1+t2)/M1
    βaprox = (b0+b1+b2+b3+b4) #β(ε=1)
    βsol1 = rr_flash_refine(K,z,βaprox)
    βsol = βsol1
    #converged to -2.2e-16 or lower, practically 0
    if abs(βsol1) < eps(_1)
        βsol =_0
    end
    if 0 <= βsol1 <= 1 #converged as expected inside range [0,1]
        βsol = βsol1
    elseif iszero(near_mean) #one βi is exacly zero or exactly one
        (rr_flash_eval(K,z,_1) < eps(_1)) && (βsol = _1)
        (rr_flash_eval(K,z,_0) < eps(_1)) && (βsol = _0)
    elseif (βsol1 < β_near0) | (βsol1 > β_near1)     #in this case, we have two asymptotes really, really near one and zero
        βaprox_near0 = near_mean
        βaprox_near1 = _1 - near_mean
        βsol_near0 = rr_flash_refine(K,z,βaprox_near0)
        βsol_near1 = rr_flash_refine(K,z,βaprox_near1)
        βsol_near0,βsol_near1
        if 0 <= βsol_near0 <= 1
            βsol = βsol_near0
        elseif 0 <= βsol_near1 <= 1
            βsol = βsol_near1
        end
    else
        βsol = _0/_0
    end
    return βsol
end
function rr_find_strongest(K,z)
    _0 = zero(first(z)+first(K))
    _1 = one(_0)
    (kmin,idmin) = findmin(K)
    (kmax,idmax) = findmax(K)
    βmin = _1/(_1-kmax)
    βmax = _1/(_1-kmin)
    idx = 0
    for i in 1:length(z)
        idx +=1
        Ki = K[i]
        idx in (idmin,idmax) && continue
        βi = _1/(_1-Ki)
        if βmin <= βi <= 0
            idmin = idx
            βmin = βi
        elseif 1 <= βi <=βmax
            idmax = idx
            βmax = βi
        end
    end

    #(4, 2, 1, 5)
    id1 = 0
    id2 = 0

    idx = 0
    strong_zβ = (-_1,-_1)
    #find the strongest components, by zi/|βi|, ignoring βmin and βmax
    for i in 1:length(z)
        Ki = K[i]
        zi = z[i]
        idx+=1
        idx in (idmin,idmax) && continue
        zβi =abs(zi*(_1-Ki))
        zβ1,zβ2 = strong_zβ
        if zβ2 <= zβi
            strong_zβ = (zβ2,zβi)
            id1 = id2
            id2 = idx
        elseif zβ1 <= zβi
            strong_zβ = (zβi,zβ2)
            id1 = idx
        end
    end

    #indices of the strongest components
    idxs = (idmin,idmax,id1,id2)
    sumz = z[idmin] + z[idmax] +z[id1] +z[id2]
    invsumz = _1/sumz
    #normalize z
    zs = (invsumz*z[idmin],invsumz*z[idmax],invsumz*z[id1],invsumz*z[id2])
    ks = (K[idmin],K[idmax],K[id1],K[id2])
    return idxs,ks,zs
    end

    function rr_vle_vapor_fraction_exact(K,z)
    #if this function is called, then we are sure that there is a solution in the interval [0,1]
    _0 = zero(first(K)+first(z))
    n = length(z)
    _1 = one(_0)
    invsumz = 1/sum(z)
    if  n == 2
        zz1,zz2 = z
        z1,z2 = zz1*invsumz,zz2*invsumz
        k1,k2 = K
        b1,b2 = 1/(1-k1),1/(1-k2)
        a1 = (z1*b2 + z2*b1)
        a0 = z1+z2
        return a1/a0
    elseif n == 3
        #0 = a0 + a1β + a2β^2
        zz1,zz2,zz3 = z
        z1,z2,z3 = zz1*invsumz,zz2*invsumz,zz3*invsumz
        k1,k2,k3 = K
        b1,b2,b3 = 1/(1-k1),1/(1-k2),1/(1-k3)
        a2 =(z1 + z2 + z3)
        a1 = -b1*(z2 + z3) - b2*(z1 + z3) - b3*(z1 + z2)
        a0 = b1*b2*z3 + b1*b3*z2 + b2*b3*z1
        Δ = a1*a1 - 4*a0*a2
        inva2 = 1/(2*a2)
        Δ2 = sqrt(Δ)*inva2
        x = -a1*inva2
        β1 = x-Δ2
        β2 = x+Δ2
        βmax = max(β1,β2)
        βmin = min(β1,β2)
        if 0 < β1 < 1
            return β1
        elseif 0 < β2 < 1
            return β2
        elseif isfinite(β1+β2)
            kmin = min(k1,k2,k3)
            kmax = max(k1,k2,k3)
            kmin > 1 && return βmax
            kmax < 1 && return βmin
        else
            return _0/_0
        end
    elseif n == 4
        #0 = a0 + a1β + a2β^2 + a3β^3
        zz1,zz2,zz3,zz4 = z
        z1,z2,z3,z4 = zz1*invsumz,zz2*invsumz,zz3*invsumz,zz4*invsumz

        k1,k2,k3,k4 = K
        b1,b2,b3,b4 = 1/(1-k1),1/(1-k2),1/(1-k3),1/(1-k4)
        a3 =(z1 + z2 + z3 + z4)
        a2 = -b1*(z2 + z3 + z4) - b2*(z1 + z3 + z4) - b3*(z1 + z2 + z4)
        a1 = b1*b2*(z3+z4) + b1*b3*(z2+z4) + b2*b3*(z1+z4) + b1*b4*(z2+z3) + b2*b4*(z1+z3) + b3*b4*(z1+z2)
        a0 = -b1*b2*b3*z4 - b1*b2*b4*z3 - b1*b3*b4*z2 - b2*b3*b4*z1
        res1 = Solvers.roots3(a0,a1,a2,a3)
        βmax = max(b1,b2,b3,b4)
        βmin = min(b1,b2,b3,b4)
        clx1,clx2,clx3 = res1
        r1,r2,r3 = real(clx1),real(clx2),real(clx3)
        rsum = r1+r2+r3
        if (r1 ≈ r2) && r1 > 1
            return r3
        elseif (r1 ≈ r3) && r1 > 1
            return r2
        elseif (r2 ≈ r3) && r2 > 1
            return r1
        end
        rmax,rmin = extrema((r1,r2,r3))
        rmid = rsum - rmax - rmin
        return rmid
    else
        return _0/_0
    end
end

function rr_flash_eval(K,z,β,non_inx = FillArrays.Fill(false,length(z)),non_iny = FillArrays.Fill(false,length(z)))
    _0 = zero(Base.promote_eltype(K,z,β))
    _1 = one(_0)
    sumz = sum(z)
    invsumz = _1/sumz
    _0βy = - 1.0 / (1.0 - β)
    _0βx = 1.0 / β
    res = _0
    res_KD2 = _0
    for i in 1:length(z)
        Ki = K[i]
        #we separate (K-1)/(1-βK) into K/(1-βK) and 1/(1-βK) for numerical reasons.
        #see;
        #Ks_eps_0 = [1.2566703532018493e-21, 3.3506275205339295, 1.0300675710905643e-23, 1.706258568414198e-39, 1.6382855298440747e-20]
        #zs_eps_0 = [0.13754371891028325, 0.29845155687154623, 0.2546683930289046, 0.08177453852283137, 0.22756179266643456]
        if Ki < eps(one(Ki))
            detKi0 = 1 - β
            detKi0 += β*Ki
            detKi = 1/detKi0
        else
            Kim1 = Ki - 1
            detKi0 = 1 + β*Kim1
            detKi = 1/detKi0
        end
        KD1 = Ki*detKi
        KD2 = -detKi
        KD = KD1 + KD2
        # modification for non-in-y components Ki -> 0
        if non_iny[i] || iszero(Ki)
            KD = _0βy
            KD2 = -_1
            KD1 = _0βy - KD2
        end
        # modification for non-in-x components Ki -> ∞
        if non_inx[i] || isinf(Ki)
            KD = _0βx
            KD1 = _0βx
            KD2 = _0
        end
        zi = z[i]
        res += zi*KD1
        res_KD2 += zi*KD2
    end
    res += res_KD2
    res *= invsumz
    return res
end

"""
    rr_flash_vapor(k,z,β)
Returns the gas phase composition, given k-values `k`, the initial molar composition `z` and the molar vapor fraction ´β´.

Each gas phase composition is calculated acording to:

    xvᵢ = kᵢzᵢ/(1+ β(kᵢ-1))
"""
function rr_flash_vapor(k,z,β)
    rr_flash_vapor!(similar(k),k,z,β)
end

function rr_flash_vapor!(y,k,z,β)
    function f(ki,zi)
    _1 = one(ki)
        return ki*zi/(_1+β*(ki-_1))
    end
    y .= f.(k,z)
    y
end

"""
    rr_flash_liquid(k,z,β)
Returns the liquid phase composition, given k-values `k`, the initial molar composition `z` and the molar vapor fraction ´β´.

Each gas phase composition is calculated acording to:

    xlᵢ = zᵢ/(1+ β(kᵢ-1))
"""
function rr_flash_liquid(k,z,β)
    rr_flash_liquid!(similar(k),k,z,β)
end

function rr_flash_liquid!(x,k,z,β)
    function f(ki,zi)
        _1 = one(ki)
        return zi/(_1+β*(ki-_1))
    end
    x .= f.(k,z)
    x
end

#obtains the minimum and maximum permissible value of β
function rr_βminmax(K,z,non_inx=FillArrays.Fill(false,length(z)), non_iny=FillArrays.Fill(false,length(z)))
    _1 = one(eltype(K))*one(eltype(z))
    βmin = Inf*_1
    βmax = -Inf*_1
    _0 = zero(_1)
    #βmin = max(0., minimum(((K.*z .- 1) ./ (K .-  1.))[K .> 1]))
    #βmax = min(1., maximum(((1 .- z) ./ (1. .- K))[K .< 1]))
    sumz = sum(z)
    for i in eachindex(K)
        Ki,zi = K[i],z[i]/sumz
        if non_inx[i] #noncondensables
            Ki = Inf*one(Ki)
        end

        if non_iny[i] # #nonvolatiles
            Ki = zero(Ki)
        end
        if Ki > 1
            # modification for non-in-x components Ki -> ∞
            not_xi = non_inx[i] || isinf(Ki)
            βmin_i = not_xi ? one(Ki)*zi : (Ki*zi - 1)/(Ki - 1)
            βmin = min(βmin,βmin_i)
        end
        if Ki < 1
            # modification for non-in-y components Ki -> 0
            not_yi = non_iny[i] || iszero(Ki)
            βmax_i = not_yi ? (1 - zi) : (zi - 1)/(Ki - 1)
            βmax = max(βmax,βmax_i)
        end
    end
    βmin = max(βmin,_0) #comment to enable negative flashes
    βmax = min(βmax,_1) #comment to enable negative flashes
    return βmin,βmax
end

function rachfordrice_β0(K,z,β0 = nothing,non_inx=FillArrays.Fill(false,length(z)), non_iny=FillArrays.Fill(false,length(z)))
    g0 = Clapeyron.rr_flash_eval(K,z,0,non_inx,non_iny)
    g1 = Clapeyron.rr_flash_eval(K,z,1,non_inx,non_iny)
    isnan(g0) && return g0,false,(g0,g0),(g0,g1)
    singlephase = false
    _1 = one(g1)
    _0 = zero(g1)

    if g0 < 0
        β = _0 #comment to enable negative flashes
        singlephase = true
    elseif g1 > 0
        β = _1 #comment to enable negative flashes
        singlephase = true
    end
    βmin,βmax = rr_βminmax(K,z,non_inx,non_iny)
    if singlephase
        βmax = max(zero(βmax),βmax)
        βmin = min(βmin,oneunit(βmin))
    end

    if β0 !== nothing
        β = β0
    else
        β = (βmax + βmin)/2
    end
    return β,singlephase,(βmin,βmax),(g0,g1)
end


#refines a rachford-rice result via Halley iterations
function rr_flash_refine(K,z,β0,non_inx=FillArrays.Fill(false,length(z)), non_iny=non_inx,limits = rr_βminmax(K,z,non_inx,non_iny))
    βmin,βmax = limits

    if βmin < 0 < 1 < βmax
        βmin = zero(βmin)
        βmax = oneunit(βmax)
    end

    _0 = zero(first(z)+first(K)+first(β0))
    _1 = one(_0)
    sumz = sum(z)
    invsumz = _1/sumz
    β = β0
    function FO(β̄)
        _0βy = - 1. / (1. - β̄)
        _0βx = 1. / β̄
        res,∂res,∂2res = _0,_0,_0
        res_KD2 = _0
        for i in 1:length(z)
            Ki = K[i]
            #we separate (K-1)/(1-βK) into K/(1-βK) and 1/(1-βK) for numerical reasons.
            #see;
            #Ks_eps_0 = [1.2566703532018493e-21, 3.3506275205339295, 1.0300675710905643e-23, 1.706258568414198e-39, 1.6382855298440747e-20]
            #zs_eps_0 = [0.13754371891028325, 0.29845155687154623, 0.2546683930289046, 0.08177453852283137, 0.22756179266643456]
            if Ki < eps(one(Ki))
                detKi0 = 1 - β̄
                detKi0 += β̄*Ki
                detKi = 1/detKi0
            else
                Kim1 = Ki - 1
                detKi0 = 1 + β̄*Kim1
                detKi = 1/detKi0
            end
            KD1 = Ki*detKi
            KD2 = -detKi
            KD = KD1 + KD2
            # modification for non-in-y components Ki -> 0
            if non_iny[i] || iszero(Ki)
                KD = _0βy
                KD2 = -_1
                KD1 = _0βy - KD2
            end
            # modification for non-in-x components Ki -> ∞
            if non_inx[i]|| isinf(Ki)
                KD = _0βx
                KD1 = _0βx
                KD2 = _0
            end
            zi = z[i]
            res += zi*KD1
            res_KD2 += zi*KD2
            ∂res -= zi*KD^2
            ∂2res += 2*zi*KD^3
        end
        res += res_KD2
        res *= invsumz
        ∂res *= invsumz
        ∂2res *= invsumz

        return res,res/∂res,∂res/∂2res
    end
    if !isfinite(β)
        return β
    end
    prob = Roots.ZeroProblem(FO,(βmin,βmax,β))
    return Roots.solve(prob,Roots.BracketedHalley())
end

function material_balance_rr_converged(w,z,β::Number,n = sum(z),ztol = sqrt(eps(eltype(β))))
    βx = SVector(1 - β,β)
    return material_balance_rr_converged(w,z,βx,n,ztol)
end

function material_balance_rr_converged(w,z,β::AbstractVector,n = sum(z),ztol = sqrt(eps(eltype(β))))
    rk = zero(Base.promote_eltype(w[1],z,β))
    np = length(β)
    #material balance test from https://github.com/WhitsonAS/Rachford-Rice-Contest
    for i in 1:length(z)
        rk1 = zero(rk)
        rk2 = zero(rk)
        zi = z[i]/n
        for j in 1:np
            Vi = β[j]*w[j][i]
            rk1 +=  Vi
            rk2 += abs(Vi)
        end
        rk = max(rk,abs(rk1 - zi)/abs(rk2 + zi))
    end
    return rk <= ztol
end


#=
it = 0
error_β = _1
error_FO = _1
while error_β > 1e-8 && error_FO > 1e-8 && it < 30
    it = it + 1
    _0βy = - 1. / (1. - β)
    _0βx = 1. / β
    res,∂res,∂2res = _0,_0,_0
    FO,dFO,d2FO = _0,_0,_0
    for i in 1:length(z)
        Kim1 = K[i] - _1
        KD = Kim1/(1+β*Kim1)
        # modification for non-in-y components Ki -> 0
        if non_iny[i]
            KD = _0βy
        end
        # modification for non-in-x components Ki -> ∞
        if non_inx[i]
            KD = _0βx
        end
        FO_i = KD
        zFOi = z[i]*FO_i
        zFOi2 = zFOi*FO_i
        zFOi3 = zFOi2*FO_i
        FO += zFOi
        dFO -= zFOi2
        d2FO += 2*zFOi3
    end
    dβ = - (2*FO*dFO)/(2*dFO^2-FO*d2FO)
    # restricted β space
    if FO < 0.
        βmax = β
    elseif FO > 0.
        βmin = β
    end

    #updatind β
    βnew = β + dβ
    if βmin < βnew && βnew < βmax
        β = βnew
    else
        dβ = (βmin + βmax) / 2 - β
        β = dβ + β
    end
    error_β = abs(dβ)
    error_FO = abs(FO)
end

return β =#