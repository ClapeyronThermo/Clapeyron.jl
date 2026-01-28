"""
    det_22(a,b,c,d)

Calculates `a*b - c*d` with less rounding error than doing it naively.
"""
function det_22(a,b,c,d)
    t = c*d
    e = fma(c,d,-t) #cd - cd
    f = fma(a,b,-t) #ab - cd
    return f-e  #ab - cd + cd - cd
end

function solve_cubic_eq(poly::AbstractVector{T}) where {T<:Real}
    tup = (poly[1],poly[2],poly[3],poly[4])
    return solve_cubic_eq(tup)
end

function depress_cubic(coeffs::NTuple{4,T}) where T
    third = one(T)/3
    a1  = one(T) / coeffs[4]
    E1,E2,E3  = -coeffs[3]*a1, coeffs[2]*a1, -coeffs[1]*a1
    E12 = E1*E1
    q_poly = (E3,-third*E2,zero(T),T(2/27))
    q = -evalpoly(E1,q_poly)
    if true || abs(q) < 10*eps(typeof(q))
        q = -det_22(2*E1,E12,9*E1,E2)/27 - E3
    end
    p = -det_22(E1,E1,3,E2)/3
    return p,q
end

function solve_cubic_eq(poly::NTuple{4,T}) where {T<:Real}
    # copied from PolynomialRoots.jl, adapted to be AD friendly
    # Cubic equation solver for complex polynomial (degree=3)
    # http://en.wikipedia.org/wiki/Cubic_function   Lagrange's method
    # poly = (a,b,c,d) that represents a + bx + cx3 + dx4

    _1 = one(T)
    third = _1/3
    a1  = one(T) / poly[4]
    E1,E2,E3  = -poly[3]*a1, poly[2]*a1, -poly[1]*a1
    s0 = E1
    p,q = depress_cubic(poly)
    A,B = -27q,-3p
    # quadratic equation: z^2 - Az + B^3=0  where roots are equal to s1^3 and s2^3
    Δ2p = (4*E2*E2*E2 + 27*E3*E3,-18*E2*E3,-E2*E2,4*E3)
    Δ2 = 27*evalpoly(E1,Δ2p) #TODO: compensated arithmetic evalpoly really needed
    if Δ2 < 10*eps(typeof(A))
        E1p = (one(E1),E1,E1*E1,E1*E1*E1)
        Δ2 = 27*dot(E1p,Δ2p)
    end
    Δ = Base.sqrt(complex(Δ2))
    @show Δ
    #Δ = (A*A - 4*B*B*B)^0.5
    if real(A*Δ)>=0 # scalar product to decide the sign yielding bigger magnitude
        s10 = 0.5 * (A + Δ)
    else
        s10 = 0.5 * (A - Δ)
    end
    r = abs(s10)
    θ = angle(s10)
    @show θ
    @show r
    s1 =  cbrt(r) * cis(θ * third)
    if s1 == 0
        s2 = s1
    else
        s2 = B / s1
    end
    zeta1 = complex(-0.5, sqrt(T(3.0))*0.5)
    zeta2 = conj(zeta1)
    return (third*(s0 + s1 + s2), third*(s0 + s1*zeta2 + s2*zeta1), third*(s0 + s1*zeta1 + s2*zeta2))
end

function solve_real_cubic_eq(poly::NTuple{4,T}) where {T<:Real}
    # copied from PolynomialRoots.jl, adapted to be AD friendly
    # Cubic equation solver for complex polynomial (degree=3)
    # http://en.wikipedia.org/wiki/Cubic_function   Lagrange's method
    # poly = (a,b,c,d) that represents a + bx + cx3 + dx4

    _1 = one(T)
    third = _1/3
    a1  = one(T) / poly[4]
    E1,E2,E3  = -poly[3]*a1, poly[2]*a1, -poly[1]*a1
    s0 = E1
    p,q = depress_cubic(poly)
    A,B = -27q,-3p
    @show A,B
    # quadratic equation: z^2 - Az + B^3=0  where roots are equal to s1^3 and s2^3
    #Δ2 = A*A - 4*B*B*B
    Δ2p = (4*E2*E2*E2 + 27*E3*E3,-18*E2*E3,-E2*E2,4*E3)
    Δ2 = 27*evalpoly(E1,Δ2p)
    if Δ2 < 10*eps(typeof(A))
        E1p = (one(E1),E1,E1*E1,E1*E1*E1)
        Δ2 = 27*dot(E1p,Δ2p)
    end
    @show Δ2
    @show A*A - 4*B*B*B
    if Δ2 > 0 #1 root only
        Δr = sqrt(Δ2)
        if A >= 0
            s10 = 0.5 * (A + Δr)
        else
            s10 = 0.5 * (A - Δr)
        end
        s1 = cbrt(s10)
        if iszero(primalval(s1))
            s2 = s1
        else
            s2 = B / s1
        end
        z1 = third*(s0 + s1 + s2)
        return z1,z1,z1
    end

    #2 or 3 roots
    
    Δ = Base.sqrt(complex(Δ2))
  
    #Δ = (A*A - 4*B*B*B)^0.5
    if real(A*Δ)>=0 # scalar product to decide the sign yielding bigger magnitude
        s10 = 0.5 * (A + Δ)
    else
        s10 = 0.5 * (A - Δ)
    end
    r = abs(s10)
    θ = angle(s10)
    s1 = cbrt(r) * cis(θ * third)
    if s1 == 0
        s2 = s1
    else
        s2 = B / s1
    end
    zeta1 = complex(-0.5, sqrt(T(3.0))*0.5)
    zeta2 = conj(zeta1)
    k1 = s1*zeta2 + s2*zeta1
    k2 = s1*zeta1 + s2*zeta2
    c1 = real(third*(s0 + s1 + s2))
    c2 = real(third*(s0 + s1*zeta2 + s2*zeta1))
    c3 = real(third*(s0 + s1*zeta1 + s2*zeta2))
    Zmin = min(c1,c2,c3)
    Zmax = max(c1,c2,c3)
    Zmid = max(min(c1,c2),min(max(c1,c2),c3))
    #refinement (https://sci-hub.st/10.1021/ie2023004)
    #=
    @show Zmin,Zmid,Zmax
    a  = poly[1]*a1
    b  = poly[2]*a1
    c  = poly[3]*a1

    for i in 1:0
        Zmin_old,Zmid_old,Zmax_old = Zmin,Zmid,Zmax
        Zmax = 0.5*Zmax -0.5*a/(Zmin*Zmid)
        Zmid = 0.5*Zmid + 0.5*(b - Zmax*Zmin)/(Zmax + Zmin)
        Zmin = 0.5*Zmin - 0.5*c - 0.5*(Zmax + Zmid)
        if abs(1 - Zmin_old/Zmin) < 8eps(typeof(Zmin))
            break
        end
    end =#
    return Zmin,Zmid,Zmax

end


"""
    roots3(pol)

Solves a cubic equation of the form pol[1] + pol[2]*x + pol[3]*x^2 + pol[4]*x^3
"""
function roots3(pol)
    a,b,c,d = pol[1],pol[2],pol[3],pol[4]
    return SVector(solve_cubic_eq((a,b,c,d)))
end

function roots3(a,b,c,d)
    x = (a,b,c,d)
    return roots3(x)
end

"""
    real_roots3(pol::NTuple{4,T}) where {T<:Real}

Given a cubic real polynom of the form `pol[1] + pol[2]*x + pol[3]*x^2 + pol[4]*x^3`,
returns `(n, zl, zg)` where `n` is the number of real roots and:
- if `n == 1`, `zl` and `zg` are equal to the only real root, the other two are complex.
- if `n == 2`, `zl` is the single real root and `zg` is the double (degenerate) real root.
- if `n == 3`, `zl` is the lowest real root and `zg` the greatest real root.

!!! info
    If there is a single root triply degenerate, e.g. with `pol == (1,3,3,1)` corresponding
    to `(x+1)^3`, this will return `n == 2` and `zl == zg` equal to the root.
"""
function real_roots3(pol::NTuple{4,T}) where {T<:Real}
    x1, x2, x3 = solve_real_cubic_eq(pol)
    if x1 == x2 == x3
        return 1,x1,x2,x3
    elseif x1 == x2 || x2 == x3
        xmin,xmax = minmax(x1,x3)
        return 2,xmin,xmax,xmax
    elseif x1 == x3
        xmin,xmax = minmax(x1,x2)
        return 2,xmin,xmax,xmax
    end
    mid12 = (x1 + x2)/2
    mid23 = (x2 + x3)/2
    between1 = evalpoly(mid12, pol)
    between2 = evalpoly(mid23, pol)
    if abs(between1) <  eps(typeof(between1)) && isapprox(x1,x2)
        (2, x3, mid12, mid12) # first the single root, then the double root
    elseif abs(between2) < eps(typeof(between2)) && isapprox(x2,x3)
        (2, x1, mid23, mid23) # first the single root, then the double root
    else
        sign1 = signbit(between1)
        sign2 = signbit(between2)
        if sign1 == sign2 # only one root
            if sign1 ⊻ (pol[4] > 0)
                (1, x1, x1, x1)
            else
                (1, x3, x3, x3)
            end
        else # three distinct roots
            (3, x1, x2, x3)
        end
    end
end

real_roots3(a,b,c,d) = real_roots3((a,b,c,d))

"""
    polyder(poly)

Returns the coefficients of the derivative of the polynomial.

"""
function polyder(x::NTuple{N,T}) where {N,T}
    return ntuple(i->x[i+1]*i,Val{N-1}())
end

#we suppose that there is a translation: xx = x + x0
function hermite5_poly(x0,x1,f0,f1,df0,df1,d2f0,d2f1)
    #Δx0 = x - x0
    #Δx1 = x - x1
    Δx10 = (x1 - x0)
    #Δx03 = Δx0^3

    Δx102 = Δx10^2
    divx = 1/Δx10
    p0 = f0
    p1 = df0

    p2 = (1//2)*d2f0
    p3 = (f1 - f0 - df0*Δx10 - (1//2)*d2f0*Δx102)*divx^3

    z4 = (3*f0 - 3*f1 + 2*(df0 + (1//2)*df1)*Δx10 + (1//2)*d2f0*Δx102)*divx^4 # * Δx03 * Δx1
    #x^3 * (x -Δx10) = x^4 - -Δx10*x^3
    p3 += -Δx10*z4
    p4 = z4

    z5 = (6*f1 - 6*f0 - 3*(df0 + df1)*Δx10 + (1//2)*(d2f1 - d2f0)*Δx102)*divx^5# * Δx03 * Δx1
    #x^3 * (x -Δx10)^2 = x^5 -2*Δx10*x^4 +Δx10*Δx10*x^3
    p3 += Δx10^2*z5
    p4 += -2*Δx10*z5
    p5 = z5
    return p0,p1,p2,p3,p4,p5
end
#0.00012669209195135698, 9.69556e-5 + 0.00012669209195135698
"""
    hermite5_poly(f,x0,x1)
    hermite5_poly(x0,x1,f0,f1,df0,df1,d2f0,d2f1)

Returns a quintic hermite polynomial, that interpolates `f` between `x0` and `x1`, using first and second derivative information.
The polynomial is translated, so that the zero is at `x0`.

## Example
```
f(x) = exp(0.2 + 3/x)
x0,x1 = 0.2,0.3
poly = hermite3_poly(f,x0,x1)
evalpoly(0.0,poly) ≈ f(x0) #true
evalpoly(x1 - x0,poly) ≈ f(x1) #true
```
"""
function hermite5_poly(f,x0,x1)
    f0,df0,d2f0 = f∂f∂2f(f,x0)
    f1,df1,d2f1 = f∂f∂2f(f,x1)
    return hermite5_poly(x0,x1,f0,f1,df0,df1,d2f0,d2f1)
end


"""
    hermite3_poly(f,x0,x1)
    hermite3_poly(x0,x1,f0,f1,df0,df1)

Returns a cubic hermite polynomial, that interpolates `f` between `x0` and `x1`, using first derivative information.
The polynomial is translated, so that the zero is at `x0`.

## Example
```
f(x) = exp(0.2 + 3/x)
x0,x1 = 0.2,0.3
poly = hermite3_poly(f,x0,x1)
evalpoly(0.0,poly) ≈ f(x0) #true
evalpoly(x1 - x0,poly) ≈ f(x1) #true
```
"""
function hermite3_poly(x0,x1,f0,f1,df0,df1)
    #=
    (2t3 - 3t2 + 0t + 1)f0 +
    (1t3 - 2t2 + 1t + 0)df0 +
    (-2t3 +3t2 + 0t + 0)f1 +
    (1t3 - 1t2 + 0t +  0)df1
    =#
    Δx⁻¹ = 1/(x1 - x0)
    p0 = f0
    p1 = df0*Δx⁻¹
    p2 = (-3*f0 -2*df0 + 3*f1 - 1*df1)*Δx⁻¹*Δx⁻¹
    p3 = (2*f0 + df0 -2*f1 + df1)*Δx⁻¹*Δx⁻¹*Δx⁻¹
    return (p0,p1,p2,p3)
end

function hermite3_poly(f,x0,x1)
    f0,df0 = f∂f(f,x0)
    f1,df1 = f∂f(f,x1)
    return hermite3_poly(x0,x1,f0,f1,df0,df1)
end

function real_roots2(coeffs)
    c,b,a = coeffs
    Δ = det_22(b,b,4*a,c)
    z0 = -b/(2*a)
    if Δ <= 0
        return 1,z0,z0
    else
        Δ2 = Base.sqrt(Δ)
        dz = Δ2/(2*a)
        if b >= 0
            x1,x2 = minmax(z0 - dz,2*c/(-b - Δ2))
            return 2,x1,x2
        else
            x1,x2 = minmax(z0 + dz,2*c/(-b + Δ2))
            return 2,x1,x2
        end
    end
end

function real_roots_qr(poly,static = false)
    p0,p1,p2,p3 = poly
    p3_inv = 1/p3
    d1 = p2*p3_inv
    d2 = p1*p3_inv
    d3 = p0*p3_inv
    M = [-d1 -d2 -d3
        1.0 0.0 0.0
        0.0 1.0 0.0]
    display(M)
    if static
        return roots3_static(M)
    else
        return eigvals(M)
    end
end

function roots3_static(A)
    #implementation of
    #NUMERICALLY STABLE EVALUATION OF CLOSED-FORM EXPRESSIONS FOR EIGENVALUES OF 3 × 3 MATRICES
    #https://arxiv.org/abs/2511.00292

    #I1 invariant
    I1 = A[1,1] + A[2,2] + A[3,3]

    #I2 invariant
    d0 = A[1,1] - A[2,2]
    d1 = A[1,1] - A[3,3]
    d2 = A[2,2] - A[3,3]
    offdiag_J2_vec = (A[1,2]*A[2,1],A[1,3]*A[3,1],A[2,3]*A[3,2])
    offdiag_J2 = sum(offdiag_J2_vec)
    diag_J2 = (d0^2 + d1^2 + d2^2)/6
    J2 = offdiag_J2 + diag_J2
    S = A - (tr(A)/3)*LinearAlgebra.I(3)
    J22 = 0.5*tr(S*S)
    @show J2-J22
    #J3 invariant
    t1 = d1 + d2
    t2 = d0 - d2
    t3 = -d0 - d1
    offdiag_J3 = A[1,2]*A[2,3]*A[3,1] + A[1,3]*A[2,1]*A[3,2]
    mixed_J3  =dot(offdiag_J2_vec,(t1,t2,t3))/3
    diag_J3 = t1*t2*t3/27
    J3 = offdiag_J3 + mixed_J3 - diag_J3
    J33 = det(S)
    #@show J3-J33
    #discriminant invariant
    #Δ = 4*J2*J2*J2 - 27*J3*J3
    #Δ = 4/27B^3 - A*A/27
    u = eig3_DX2(A)
    v = eig3_DX2(transpose(A))
    Δ = zero(J3)
    w = (9,6,6,6,8,8,8,2,2,2,2,2,2,1)

    for i in 1:14
        @show u[i]*v[i]*w[i]
        Δ += u[i]*v[i]*w[i]
    end
    #@show Δ
    #@show 4*J2^3 - 27*J3^2
    A = J3*27
    B = J2*3
    #@show A,B
    Δ = 4*J2*J2*J2 - 27*J3*J3
    #@show (A*A - 4*B*B*B)
    t = sqrt(27*complex(Δ))/(27*J3)
    ϕ = atan(t)
    k = 2*sqrt(3*complex(J2))
    λ1 = I1 + k*cos((ϕ + 2π)/3)
    λ2 = I1 + k*cos((ϕ + 4π)/3)
    λ3 = I1 + k*cos(ϕ/3)
    return [λ1/3,λ2/3,λ3/3]
end

function eig3_DX2(A)
    A00,A11,A22 = A[1,1],A[2,2],A[3,3]
    A01,A02,A12 = A[1,2],A[1,3],A[2,3]
    A10,A20,A21 = A[2,1],A[3,1],A[2,3]
    d0 = A00 - A11
    d1 = A00 - A22
    d2 = A11 - A22
    r0 = A01 * A12 * A20 - A02 * A10 * A21
    r1 = -A01 * A02 * d2 + A01 * A01 * A12 -
                A02 * A02 * A21
    r2 = A01 * A21 * d1 - A01 * A01 * A20 +
                A02 * A21 * A21
    r3 = A02 * A12 * d0 + A01 * A12 * A12 -
                A02 * A02 * A10
    r4 = A01 * A12 * d1 - A01 * A02 * A10 +
                A02 * A12 * A21
    r5 = A02 * A21 * d0 - A01 * A02 * A20 +
                A01 * A12 * A21
    r6 = -A02 * A10 * d2 + A01 * A10 * A12 -
                A02 * A12 * A20
    r7 = A12 * d0 * d1 - A02 * A10 * d1 +
                A01 * A10 * A12 - A12 * A12 * A21
    r8 = A12 * d0 * d1 - A02 * A10 * d0 +
                A02 * A12 * A20 - A12 * A12 * A21
    r9 = A01 * d1 * d2 + A02 * A21 * d2 +
                A01 * A02 * A20 - A01 * A01 * A10
    r10 = A01 * d1 * d2 + A02 * A21 * d1 +
                A01 * A12 * A21 - A01 * A01 * A10
    r11 = -A02 * d0 * d2 + A01 * A12 * d0 +
                A02 * A12 * A21 - A02 * A02 * A20
    r12 = A02 * d0 * d2 + A01 * A12 * d2 -
                A01 * A02 * A10 + A02 * A02 * A20
    r13 = d0 * d1 * d2 - A01 * A10 * d0 + A02 * A20 * d1 -
                A12 * A21 * d2
    return (r0,r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13)
end

function eig3_DX(A)
    A00,A11,A22 = A[1,1],A[2,2],A[3,3]
    A01,A02,A12 = A[1,2],A[1,3],A[2,3]
    A10,A20,A21 = A[2,1],A[3,1],A[2,3]
    d0 = A00 - A11
    d1 = A00 - A22
    d2 = A11 - A22
    r1 = A01*A12*A20 − A02*A10*A21
    r2 = −A01*A02*d2 +(A01^2)*A12 −(A02^2)*A21
    r3 = A01*A21*d1 −(A01^2)*A20 + A02*(A21^2)
    r4 = A02*A12*d0 + A01*(A12^2) −(A02^2)*A10
    r5 = A01*A12*d1 − A01*A02*A10 + A02*A12*A21
    r6 = A02*A21*d0 − A01*A02*A20 + A01*A12*A21
    r7 = −A02*A10*d2 + A01*A10*A12 − A02*A12*A20
    r8 = A12*d0*d1 − A02*A10*d1 + A01*A10*A12 −(A12^2)*A21
    r9 = A12*d0*d1 − A02*A10*d0 + A02*A12*A20 −(A12^2)*A21
    r10 = A01*d1*d2 − A02*A21*d2 + A01*A02*A20 −(A01^2)*A10
    r11 = A01*d1*d2 + A02*A21*d1 + A01*A12*A21 −(A01^2)*A10
    r12 = −A02*d0*d2 + A01*A12*d0 + A02*A12*A21 −(A02^2)*A20
    r13 = A02*d0*d2 + A01*A12*d2 − A01*A02*A10 +(A02^2)*A20
    r14 = d0*d1*d2 − A01*A10*d0 + A02*A20*d1 − A12*A21*d2
    r = (r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14)
    return r
end 