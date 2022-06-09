"""
    rr_vle_vapor_fraction(K,z,α = NaN)

Given a vector of K values and a vector of compositions, calculates the vapor fraction `β`.
the algorithm is a modification of _(1)_ , with safeguards for extreme cases.

If the algorithm fails to converge, returns `NaN`. if it converges to a value `β ∉ [0,1]`, returns `-Inf` or `Inf`, depending on the case.

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
    if n_z <= 4
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
    βsol1 =  rr_flash_refine(K,z,βaprox)
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
        βsol_near0 =  rr_flash_refine(K,z,βaprox_near0) 
        βsol_near1 =  rr_flash_refine(K,z,βaprox_near1) 
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
    (kmin,idmin) =  findmin(K)
    (kmax,idmax) =  findmax(K)
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
    if  n == 2
        z1,z2 = z
        k1,k2 = K
        b1,b2 = 1/(1-k1),1/(1-k2)
        a1 = (z1*b2 + z2*b1)
        a0 = z1+z2
        return a1/a0
    elseif n == 3
        #0 = a0 + a1β + a2β^2
        z1,z2,z3 = z
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
        z1,z2,z3,z4 = z
        k1,k2,k3,k4 = K
        b1,b2,b3,b4 = 1/(1-k1),1/(1-k2),1/(1-k3),1/(1-k4)
        a3 =(z1 + z2 + z3 + z4)
        a2 = -b1*(z2 + z3 + z4) - b2*(z1 + z3 + z4) - b3*(z1 + z2 + z4)
        a1 = b1*b2*(z3+z4) + b1*b3*(z2+z4) + b2*b3*(z1+z4) + b1*b4*(z2+z3) + b2*b4*(z1+z3) + b3*b4*(z1+z2)
        a0 = -b1*b2*b3*z4 - b1*b2*b4*z3 - b1*b3*b4*z2 - b2*b3*b4*z1
        res1 =  Solvers.roots3(a0,a1,a2,a3)
        βmax = max(b1,b2,b3,b4)
        βmin = min(b1,b2,b3,b4)
        clx1,clx2,clx3 = res1
        r1,r2,r3 = real(clx1),real(clx2),real(clx3)
        rsum = r1+r2+r3
        rmax,rmin = extrema((r1,r2,r3))
        rmid = rsum - rmax - rmin
        return rmid
    else
        return _0/_0
    end
    end
    
function rr_flash_eval(K,z,β,normalize=true)
        _0 = zero(first(z)+first(K)+first(β))
        _1 = one(_0)
        if normalize
            sumz = sum(z)
            invsumz = _1/sumz
        else
            invsumz = _1
        end
        res = _0
        for i in 1:length(z)
            Kim1 = K[i] - _1
            res += invsumz*z[i]*Kim1/(1+β*Kim1)
        end
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

function rr_flash_refine(K,z,β) 
    function FO(β̄ )
        _0 = zero(first(z)+first(K)+first(β̄ ))
        _1 = one(_0)
        sumz = sum(z)
        invsumz = _1/sumz
        res,∂res,∂2res = _0,_0,_0
        for i in 1:length(z)
            Kim1 = K[i] - _1
            KD = Kim1/(1+β̄ *Kim1)
            res += invsumz*z[i]*KD
            ∂res -= invsumz*z[i]*KD^2
            ∂2res += 2*invsumz*z[i]*KD^3
        end
        return res,∂res,∂2res
    end
    return Solvers.halley(FO,β)
end