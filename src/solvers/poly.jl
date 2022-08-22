function evalpoly_cheb(x::S, ch::AbstractVector{T}) where {T,S}
    R = promote_type(T, S)
    l = length(ch)
    l == 0 && return zero(R)
    l == 1 && return R(ch[1])
    c0 = ch[l - 1]
    c1 = ch[l]
    for i in (l-2):-1:1
        c0, c1 = ch[i] - c1, c0 + c1 * 2x
    end
    return R(c0 + c1 * x)
end


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
# poly = (a,b,c,d) that represents ax3 + bx2 + cx2 + d 

_1 = one(T)
third = _1/3
a1  =  complex(one(T) / poly[4])
E1  = -poly[3]*a1
E2  =  poly[2]*a1
E3  = -poly[1]*a1
s0  =  E1
E12 =  E1*E1
A   =  det_22(2*E1,E12,9*E1,E2) + 27*E3
#A   =  2*E1*E12 - 9*E1*E2 + 27*E3 # = s1^3 + s2^3
B = det_22(E1,E1,3,E2)
#B   =  E12 - 3*E2                 # = s1 s2
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
solves a cubic equation of the form pol[1] + pol[2]*x + pol[3]*x^2 + pol[4]*x^3
"""
function roots3(pol) 
return SVector(solve_cubic_eq(pol))
end

function roots3(a,b,c,d) 
x = (a,b,c,d)
return roots3(x)
end