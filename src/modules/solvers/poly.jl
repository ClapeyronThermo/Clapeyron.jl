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


function __roots2_det(a::T, b::T, c::T) where T
    ā = abs(a)
    b̄ = abs(b)
    c̄ = abs(c)
    sa = a > 0
    sc = c > 0
    x = sqrt(ā) * sqrt(c̄)
    _0 = zero(x)
    if sa == sc
        # Case ac > 0: need sqrt(b^2 - ac)
        # sqrt(b2 - ac) = sqrt(b - sqrt(ac))*sqrt(b + sqrt(ac))
        # Compute ac/x and its error compensation.
        if ā > x
            ax = ā / x
            ε_ax = (fma(ax, x, -ā)) / x
            acx = ax * c̄
            ε_acx = fma(ax, -c̄, acx) + ε_ax * c̄
        else
            cx = c̄ / x
            ε_cx = (fma(cx, x, -c̄)) / x
            acx = ā * cx
            ε_acx = fma(-ā, cx, acx) + ε_cx * ā
        end
        #d* = |b| - sqrt(ac) - (ac/x - x)/2 with error correction
        d̄ = b̄ - x - (acx - x - ε_acx) / 2

        if iszero(d̄)
            return true,_0
        else
            d_is_real = d̄ > 0
            return d_is_real,sqrt(abs(d̄)) * sqrt(b̄ + x)
        end
    else
        # Case ac < 0: need sqrt(b^2 + ac) = sqrt(b^2 + |ac|)
        # Use hypot‑style ordering to avoid overflow.
        if b̄ > x
            z = x / b̄
            return true,b̄ * hypot(one(_0),z)
        else
            z = b̄ / x
            return true,x * hypot(one(_0),z)
        end
    end
end

function __roots2(pol::NTuple{3,T}) where T
    c,b,a = pol
    return __roots2(a,b,c)
end

function __roots2(a::T, b::T, c::T) where T
    # Solve a*x^2 + b*x + c = 0, return real roots as a vector.
    # Corresponds to quadratic-solutions in Racket.
    
    # Handle the linear case a == 0
    if iszero(a)
        r0 = -c / b
        return true, r0,r0
    end

    # Use b/2 for simplification
    b_half = b / 2.0
    
    # Compute discriminant sqrt((b/2)^2 - a*c)
    if a isa Integer && b isa Integer && c isa Integer
        d = (b*b)//4 - a*c
        sqrt_abs_d = sqrt(abs(d))
        d_is_real = d > 0
    else
        d_is_real,sqrt_abs_d = __roots2_det(a, b_half, c)
    end
    
    # Handle double root
    if iszero(sqrt_abs_d)
        r0 = -b_half / a
        return true,r0,r0
    end

    # Handle complex roots
    if !d_is_real
        re0 = -b_half / a
        im0 = sqrt_d / abs(a)
        return false,re0,im0
    end  

    # Use a/c swapping trick to avoid cancellation
    if b < 0
        # When b < 0, -b/2 + sqrt_d is more stable
        r1 = c / (sqrt_abs_d - b_half)
        r2 = (sqrt_abs_d - b_half) / a
    else
        # When b >= 0, use equivalent forms that avoid cancellation
        r1 = (b_half + sqrt_abs_d) / (-a)
        r2 = -c / (b_half + sqrt_abs_d)
    end
    return true, r1, r2
end



#=
cbrt(a + b) + cbrt(a - b)

this is equal to cbrt(a)
=#
function stable_cbrt_sum(a, b)
    if iszero(primalval(a)) || iszero(primalval(b))
        return cbrt(a+b) + cbrt(a-b)
    end
    n = 2 * a
    c̄ = a + b
    c̲ = a - b
    t1 = cbrt(c̄ * c̄)   # (a+b)^(2/3)
    t2 = cbrt(c̲ * c̄)   # (a^2 - b^2)^(1/3), but computed as product
    t3 = cbrt(c̲ * c̲)   # (a-b)^(2/3)
    d = t1 - t2 + t3
    return n / d
end

#=base cubic root solver
returns num_real_roots,r1,r2,r3
if num_real_roots < 3 then r2 and r3 represent the real and imaginary part of the complex roots:
    c1 = Complex(r1)
    c2 = Complex(r3,r2)
    c3 = Complex(r3,-r2)

where r2 is the imaginary part and r3 is the real part.

if num_real_roots == 2, then r2 can be ignored and you can just take r1 and r3.
=#
function __roots3(pol::NTuple{4,T}) where T
    #x³ + ax² + bx + c
    pp0,pp1,pp2,pp3 = pol
    a = pp2/pp3
    b = pp1/pp3
    c = pp0/pp3
    _0 = zero(c)
    abc = (a,b,c)
    

    I1 = a
    I1 = -a
    _3J2 = a*a - 3*b # J2 = = a²/3 - b

    _27J3 = fma(-27, c, (((9 * b) - ((a + a) * a)) * a))   #J3 = ab/3 - 2a³/27 - c
    #Δ = ((((((a * a) * b) * b) - (4 * b^3)) - ((4 * c) * a^3)) - (c * (27 * c))) + (((18 * a) * b) * c) #Δ = a² b² - 4b³ - 4a³ c - 27c² + 18abc
    Δ = fma(b, 
        (b * det_22(a,a,4,b)), 
        (c * ((-4 * a * a * a) - 
        9*det_22(3,c,2*a,b)))
        )
    if !isfinite(Δ)
        r0 = oftype(Δ,NaN)
        return 0,r0,r0,r0,Δ
    end
    
    _pol = (oftype(Δ,pp0),oftype(Δ,pp1),oftype(Δ,pp2),oftype(Δ,pp3))
    if Δ >= _0
        sqrt_3J2 = sqrt(max(_0, _3J2))
        α = 2 * sqrt_3J2
        ϕ = atan(sqrt(max(_0, 27 * Δ)), _27J3)  # atan(y,x)

        λ1 = ϕ / 3
        λ2 = (ϕ + 2π) / 3
        λ3 = (ϕ + 4π) / 3
        r1 = (I1 + α * cos(λ1)) / 3
        r2 = (I1 + α * cos(λ2)) / 3
        r3 = (I1 + α * cos(λ3)) / 3
        #refine
        r1,r2,r3 = durand_kerner_refine(_pol,(r1,r2,r3))    
        #sort
        r1,r2,r3 = (r1,minmax(r2,r3)...)
        r1,r2,r3 = (minmax(r1,r2)...,r3)
        r1,r2,r3 = (r1,minmax(r2,r3)...)

        if r1 == r2 && r2 == r3
            r1 = refine_poly(_pol,r1)
            return complex_roots_from_r1(abc,r1,Δ)
        elseif r1 == r2
            return 2,r3,_0,r1,Δ
        elseif r2 == r3
            return 2,r1,_0,r3,Δ
        end

        r12 = (r1 + r2)/2
        r23 = (r2 + r3)/2
        f12 = evalpoly(r12, _pol)
        f23 = evalpoly(r23, _pol)
        if abs(f12) <  eps(typeof(f12)) && isapprox(r1,r2)
            return (2, r3, _0, r12, Δ) # first the single root, then the double root
        elseif abs(f23) < eps(typeof(f23)) && isapprox(r2,r3)
            return (2, r1, _0, r23, Δ) # first the single root, then the double root
        else
            sign1 = signbit(f12)
            sign2 = signbit(f23)
            if sign1 == sign2 # only one root
                if sign1 ⊻ (pp3 > 0)
                    r1 = refine_poly(_pol,r1)
                    return complex_roots_from_r1(abc,r1,Δ)
                else
                    r3 = refine_poly(_pol,r3)
                    return complex_roots_from_r1(abc,r3,Δ)
                end
            else # three distinct roots
                return (3, r1, r2, r3,Δ)
            end
        end
    else
        D = -Δ / 108            # D > 0
        sqrtD = sqrt(D)
        J3_half = fma(-a, b / -6, fma(-0.5, c, (((a * a) * a) * (-37 // 999))))
        r1 = stable_cbrt_sum(J3_half,sqrtD) - a/3
        r1 = refine_poly(_pol,r1)
        return complex_roots_from_r1(abc,r1,Δ)
    end
end

function complex_roots_from_r1(abc::NTuple{3,T},r1::T,Δ::T) where T
    a,b,c = abc
    A2 = a + r1
    B2 = evalpoly(r1,(b,a,one(T))) #b + a*r1 + r1^2
    disc2 = fma(r1 + a, fma(-3, r1, a), -4 * b)
    real_part = -A2 / 2
    imag_part = sqrt(max(zero(disc2), -disc2)) / 2
    return 1,r1,imag_part,real_part,Δ
end

function refine_poly(poly::NTuple{4,T},r1::T) where T
    dpoly = polyder(poly)
    x = r1
    _1 = one(r1)
    tol = T(1e-14)
    for _ in 1:100
        f = evalpoly(x,poly)
        df = evalpoly(x,dpoly)
        if abs(df) <= 10*eps(eltype(r1)) * max(_1, abs(x))
            break
        end
        dx = f / df
        x = x - dx
        if abs(dx) <= tol * max(_1, abs(x))
            break
        end
    end
    return x
end

function durand_kerner_refine(pol::NTuple{4,T},r123::NTuple{3,T}) where T
    ω = 0.5
    r = r123
    dr1 = abs(r[1] - r[2])
    dr2 = abs(r[2] - r[3])
    dr3 = abs(r[1] - r[3])
    dr = min(dr1,min(dr2,dr3))
    for iter in 1:20
        rnew = r
        @inbounds for i in 1:3
            x = r[i]
            # Evaluate p(x)
            fx = evalpoly(x,pol)

            # Compute denominator product over all j != i
            denom = one(T)
            @inbounds for j in 1:3
                if j != i
                    denom *= (r[i] - r[j])
                end
            end

            # Avoid division by zero for (near‑)multiple roots
            dx = ω * fx / denom
            if abs(dx) > 0.5*dr
                ri = x
            else
                ri = x - dx
            end  
         
            rnew = Base.setindex(rnew,ri,i)
        end
        rn1,rn2,rn3 = rnew
        rr1,rr2,rr3 = r
        max_change = max(abs(rn1-rr1),max(abs(rn2-rr2),abs(rn3-rr3)))
        r = rnew
        if max_change < 1e-12
            break
        end
    end
    return r
end

"""
    roots3(pol)

Solves a cubic equation of the form pol[1] + pol[2]*x + pol[3]*x^2 + pol[4]*x^3
"""
function roots3(pol)
    _pol = promote(pol[1],pol[2],pol[3],pol[4])
    nr,r1,r2,r3,Δ = __roots3(_pol)
    if nr == 3 || nr == 0
        return SVector(Complex(r1),Complex(r2),Complex(r3))
    end
    return SVector(Complex(r1),Complex(r3,r2),Complex(r3,-r2))
end

function roots3(a,b,c,d)
    x = promote(a,b,c,d)
    return roots3(x)
end

"""
    real_roots3(pol::NTuple{4,T}) where {T<:Real}

Given a cubic real polynom of the form `pol[1] + pol[2]*x + pol[3]*x^2 + pol[4]*x^3`,
returns `(n, r1, r2, r3)` where `n` is the number of real roots and:
- if `n == 1`, returns `(r1,r1,r1)`, only the real root.
- if `n == 2`, returns `(r1,rx,r3)`, sorted. where `rx` will be the value of the double root.
- if `n == 3`, returns `(r1,r2,r3)`, sorted.

!!! info
    If there is a single root triply degenerate, e.g. with `pol == (1,3,3,1)` corresponding
    to `(x+1)^3`, the solver may return `2,(r1,r1,r1)`, where the double and the triple root being equal.
"""
function real_roots3(pol::NTuple{4,T}) where {T<:Real}
    nr,r1,r2,r3,Δ = __roots3(pol)
    nr == 1 && (return nr,r1,r1,r1)
    if nr == 2
        if r3 < r1
            return nr,r3,r3,r1
        else
            return nr,r1,r3,r3
        end
    end
    return nr,r1,r2,r3
end

real_roots3(a,b,c,d) = real_roots3((a,b,c,d))

function real_roots2(a,b,c)
    return real_roots2(promote(c,b,a))
end

function real_roots2(pol)
    _pol = promote(pol[1],pol[2],pol[3])
    return real_roots2(_pol)
end


function real_roots2(pol::NTuple{3,T}) where T
    return __roots2(pol)
end



"""
    polyder(poly)

Returns the coefficients of the derivative of the polynomial.

"""
function polyder(x::NTuple{N,T}) where {N,T}
    return ntuple(i->x[i+1]*i,Val{N-1}())
end

function polyder(x)
    xx = @view x[2:end]
    return map(p -> p[2]*(p[1]),enumerate(xx))
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