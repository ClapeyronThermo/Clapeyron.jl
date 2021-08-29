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
    
    if length(z) <= 4
        return rr_vle_vapor_fraction_exact(K,z)
    end
    
    #z_ℵ1 is normalized
    ℵ1,k_ℵ1,z_ℵ1 = Clapeyron.rr_find_strongest(K,z)
    β0 = rr_vle_vapor_fraction_exact(k_ℵ1,z_ℵ1)
    b0 = β0
    _0 = zero(first(z)+first(K))
    _1 = one(_0)
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
    #@show b0,b1,b2,b3,b4
    βx = (b0+b1+b2+b3+b4)
    #now we have 5 points to select.
    r0 = rr_flash_eval(K,z,b0)
    β1 = b0+b1
    r1 = rr_flash_eval(K,z,β1)
    β2 = β1+b2
    r2 = rr_flash_eval(K,z,β2)
    β3 = β2+b3
    r3 = rr_flash_eval(K,z,β3)
    β4 = β3+b4
    r4 = rr_flash_eval(K,z,β4)
    bx = SA[r0,r1,r2,r3,r4]
    all_β = (β0,β1,β2,β3,β4)
    idx_b = sortperm(bx)
    sort_β = map(x->getindex(all_β,x),idx_b)
    #TODO: refine 
    return β4
end
function rr_find_strongest(K,z)
    _0 = zero(first(z)+first(K))
    _1 = one(_0)
    (_,idmin) =  findmin(K)
    (_,idmax) =  findmax(K)
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
    #=here are some tests, from the paper.
    the paper has errors, z2 has 10 components, 9 and 10 repeated.
    and after that, they dont include the bmax,bmin in their selected pair.
    z1 = [0.30,0.15, 0.05, 0.02, 0.01, 0.02, 0.02, 0.03, 0.07, 0.33]
    k1 = [3.0E+00, 2.0E+00, 1.1E+00, 8.0E-01, 4.0E-01, 1.0E-01, 5.0E-02, 2.0E-02, 1.0E-02, 1.0E-04]
    idxs1,ks1,zs1 = Clapeyron.rr_find_strongest(k1,z1)
    @test all(in.(idxs1,Ref((1,2,9,10))))
    
    z2 = [0.00034825,0.01376300,0.13084000,0.10925000,0.00001000,0.51009000,0.23564000,0.00006000]
    k2 = [5.2663E+02, 5.0400E+01, 1.6463E+00, 8.7450E-01, 1.5589E-01, 3.6588E-02, 2.6625E-02, 4.8918E-06]
    idxs2,ks2,zs2 = Clapeyron.rr_find_strongest(k2,z2)
    @test_broken all(in.(idxs2,Ref((1,2,6,7))))
    #the paper is wrong here.
    z3 = [0.0187002, 0.0243002, 0.5419054, 0.0999010, 0.0969010, 0.0400004, 0.0212002, 0.0148001, 0.0741507, 0.0350404, 0.0173602, 0.0157402]
    k3 = [1.32420, 1.12778, 1.22222, 1.11760, 9.88047E-01, 8.94344E-01, 7.87440E-01, 7.43687E-01, 8.11797E-01, 6.93279E-01, 5.09443E-01, 2.28721E-01]
    idxs3,ks3,zs3 = Clapeyron.rr_find_strongest(k3,z3)
    @test all(in.(idxs3,Ref((1,3,9,12))))
    
    =#
    function rr_vle_vapor_fraction_exact(K,z)
    n = length(z)
    if n == 2
        #0 = a1b -a0
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
        bmin = min(b1,b2,b3)
        bmax = max(b1,b2,b3)
        a2 =(z1 + z2 + z3) 
        a1 = -b1*(z2 + z3) - b2*(z1 + z3) - b3*(z1 + z2)
        a0 = b1*b2*z3 + b1*b3*z2 + b2*b3*z1
        Δ = a1*a1 - 4*a0*a2
        inva2 = 1/(2*a2)
        Δ2 = sqrt(Δ)*inva2
        x = -a1*inva2
        β1 = x-Δ2
        β2 = x+Δ2
        #@show β1,β2
        #@show bmin,bmax
        if bmin < β1 < bmax
            return β1
        else
            return β2
        end
    elseif n == 4
        #0 = a0 + a1β + a2β^2 + a3β^3
        z1,z2,z3,z4 = z
        k1,k2,k3,k4 = K
        b1,b2,b3,b4 = 1/(1-k1),1/(1-k2),1/(1-k3),1/(1-k4)
        bmin = min(b1,b2,b3,b4)
        bmax = max(b1,b2,b3,b4)
        a3 =(z1 + z2 + z3 + z4)
        a2 = -b1*(z2 + z3 + z4) - b2*(z1 + z3 + z4) - b3*(z1 + z2 + z4)
        a1 = b1*b2*(z3+z4) + b1*b3*(z2+z4) + b2*b3*(z1+z4) + b1*b4*(z2+z3) + b2*b4*(z1+z3) + b3*b4*(z1+z2)
        a0 = -b1*b2*b3*z4 - b1*b2*b4*z3 - b1*b3*b4*z2 - b2*b3*b4*z1
        res1 =  Solvers.roots3(a0,a1,a2,a3)
        clx1,clx2,clx3 = res1
        r1,r2,r3 = real(clx1),real(clx2),real(clx3)
        rsum = r1+r2+r3
        rmax,rmin = extrema((r1,r2,r3))
        rmid = rsum - rmax - rmin
        return rmid
    else
        _0 = zero(first(K)+first(z))
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
        #@show β,res

        return res
    end
     #porting those functions from my ex pkg.
   
"""
    rr_flash_vapor(k,z,β) 
Returns the gas phase composition, given k-values `k`, the initial molar composition `z` and the molar vapor fraction ´β´.

Each gas phase composition is calculated acording to:

    xvᵢ = kᵢzᵢ/(1+ β(kᵢ-1))
"""
function rr_flash_vapor(k,z,β)
    function f(ki,zi)
    _1 = one(ki) 
        return ki*zi/(_1+β*(ki-_1))
    end

    return map(f,k,z)
end

"""
    rr_flash_liquid(k,z,β) 
Returns the liquid phase composition, given k-values `k`, the initial molar composition `z` and the molar vapor fraction ´β´.

Each gas phase composition is calculated acording to:

    xlᵢ = zᵢ/(1+ β(kᵢ-1))
"""
function rr_flash_liquid(k,z,β)
    function f(ki,zi)
        _1 = one(ki) 
        return zi/(_1+β*(ki-_1))
    end
    return map(f,k,z)
end
