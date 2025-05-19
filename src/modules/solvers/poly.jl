"""
    det_22(a,b,c,d)

Calculates `a*b - c*d` with less rounding error than doing it naively
"""
function det_22(a,b,c,d)
    t = c*d
    e = muladd(c,d,-t) #cd - cd
    f = muladd(a,b,-t) #ab - cd
    return f-e  #ab - cd + cd - cd
end

function solve_cubic_eq(poly::AbstractVector{T}) where {T<:Real}
    tup = (poly[1],poly[2],poly[3],poly[4])
    return solve_cubic_eq(tup)
end

function solve_cubic_eq(poly::NTuple{4,T}) where {T<:Real}
    # copied from PolynomialRoots.jl, adapted to be AD friendly
    # Cubic equation solver for complex polynomial (degree=3)
    # http://en.wikipedia.org/wiki/Cubic_function   Lagrange's method
    # poly = (a,b,c,d) that represents a + bx + cx3 + dx4

    _1 = one(T)
    third = _1/3
    a1  = complex(one(T) / poly[4])
    E1  = -poly[3]*a1
    E2  = poly[2]*a1
    E3  = -poly[1]*a1
    s0  = E1
    E12 = E1*E1
    A   = det_22(2*E1,E12,9*E1,E2) + 27*E3
    #A   = 2*E1*E12 - 9*E1*E2 + 27*E3 # = s1^3 + s2^3
    B = det_22(E1,E1,3,E2)
    #B   = E12 - 3*E2                 # = s1 s2
    # quadratic equation: z^2 - Az + B^3=0  where roots are equal to s1^3 and s2^3
    Δ2 = det_22(A,A,4*B*B,B)
    Δ = Base.sqrt(Δ2)
    #Δ = (A*A - 4*B*B*B)^0.5
    if real(conj(A)*Δ)>=0 # scalar product to decide the sign yielding bigger magnitude
        s1 = exp(log(0.5 * (A + Δ)) * third)
    else
        s1 = exp(log(0.5 * (A - Δ)) * third)
    end
    if s1 == 0
        s2 = s1
    else
        s2 = B / s1
    end
    zeta1 = complex(-0.5, sqrt(T(3.0))*0.5)
    zeta2 = conj(zeta1)
    return (third*(s0 + s1 + s2), third*(s0 + s1*zeta2 + s2*zeta1), third*(s0 + s1*zeta1 + s2*zeta2))
end

"""
    roots3(pol)

Solves a cubic equation of the form pol[1] + pol[2]*x + pol[3]*x^2 + pol[4]*x^3
"""
function roots3(pol)
    a,b,c,d = pol[1],pol[2],pol[3],pol[4]
    return solve_cubic_eq((a,b,c,d))
end

function roots3(a,b,c,d)
    x = (a,b,c,d)
    return roots3(x)
end

"""
    real_roots3(pol::NTuple{4,T}) where {T<:Real}

Given a cubic real polynom of the form `pol[1] + pol[2]*x + pol[3]*x^2 + pol[4]*x^3`,
return `(n, zl, zg)` where `n` is the number of real roots and:
- if `n == 1`, `zl` and `zg` are equal to the only real root, the other two are complex.
- if `n == 2`, `zl` is the single real root and `zg` is the double (degenerate) real root.
- if `n == 3`, `zl` is the lowest real root and `zg` the greatest real root.

!!! info
    If there is a single root triply degenerate, e.g. with `pol == (1,3,3,1)` corresponding
    to `(x+1)^3`, this will return `n == 2` and `zl == zg` equal to the root.
"""
function real_roots3(pol::NTuple{4,T}) where {T<:Real}
    z1, z2, z3 = sort(SVector(roots3(pol)); by=real)
    x1, x2, x3 = real(z1), real(z2), real(z3)
    mid12 = (x1 + x2)/2
    mid23 = (x2 + x3)/2
    between1 = evalpoly(mid12, pol)
    between2 = evalpoly(mid23, pol)
    if abs(between1) <  eps(typeof(between1))
        (2, x3, mid12, mid12) # first the single root, then the double root
    elseif abs(between2) < eps(typeof(between2))
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

returns the coefficients of the derivative of the polynomial.

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
