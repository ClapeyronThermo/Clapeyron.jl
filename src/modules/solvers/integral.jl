const gk21_x = [-0.9956571630258080807355,-0.973906528517171720078,-0.9301574913557082260012,-0.8650633666889845107321,-0.7808177265864168970637,-0.6794095682990244062343,-0.562757134668604683339,
        -0.4333953941292471907993,-0.294392862701460198131,-0.1488743389816312108848,0,0.1488743389816312108848,0.2943928627014601981311,0.4333953941292471907993,
        0.562757134668604683339,0.6794095682990244062343,0.7808177265864168970637,0.865063366688984510732,0.9301574913557082260012,0.973906528517171720078,0.9956571630258080807355]

const gk21_w = [0.0116946388673718742781,0.0325581623079647274788,0.0547558965743519960314,0.075039674810919952767,0.093125454583697605535,0.1093871588022976418992,0.123491976262065851078,
        0.134709217311473325928,0.142775938577060080797,0.1477391049013384913748,0.149445554002916905665,0.1477391049013384913748,0.1427759385770600807971,0.134709217311473325928,
        0.123491976262065851078,0.109387158802297641899,0.093125454583697605535,0.075039674810919952767,0.05475589657435199603138,0.032558162307964727479,0.0116946388673718742781]


"""
    integral21(f,a,b)

Performs integration via a 21-point Gauss-Kronrod rule. ForwardDiff-friendly.

## example
```
julia> integral21(exp,1,1.5) - (exp(1.5) - exp(1))
0.0

julia> integral21(abs2,2,5) - (125/3 - 8/3)
-7.105427357601002e-15
```
"""
function integral21(Base.@specialize(f),a,b)
    w = gk21_w
    x = gk21_x
    res = zero(a+b)
    for i in 1:21
        xs = 0.5*(x[i] + 1)*(b - a) + a
        res += w[i]*f(xs)
    end
    return 0.5*(b - a)*res
end

# u: the quadrature points in the domain [0, ∞) for Gauss-Laguerre integration with 5 points 
# w: the weights associated with these quadrature points
# Note: Additional quadrature points can be generated with the python code at https://en.wikipedia.org/wiki/Gauss%E2%80%93Laguerre_quadrature
const laguerre5_u = (0.26356031971814109102031,1.41340305910651679221800,3.59642577104072208122300,7.08581000585883755692200,12.6408008442757826594300)
const laguerre5_w = (0.5217556105828086524759,0.3986668110831759274500,7.5942449681707595390e-2,3.6117586799220484545e-3,2.3369972385776227891e-5)

const laguerre10_u = (0.13779347054049243, 0.7294545495031705, 1.808342901740316, 3.4014336978548996, 5.552496140063804, 8.330152746764497, 11.843785837900066, 16.279257831378104, 21.99658581198076, 29.92069701227389)
const laguerre10_w = (0.30844111576502015, 0.40111992915527356, 0.2180682876118094, 0.062087456098677746, 0.0095015169751811, 0.0007530083885875388, 2.8259233495995656e-5, 4.2493139849626863e-7, 1.8395648239796308e-9, 9.911827219609008e-13)

"""
    laguerre5(f,r = 1, a = 0)
    
Performs a 5-point translated gauss-laguerre integration of the form:

```
∫exp(-ry)f(y) dy,    y ∈ (a,∞)
```
"""
function laguerre5(Base.@specialize(f),r = 1.,a = 0.)
   return _laguerrex(f,r,a,laguerre5_u,laguerre5_w)
end

"""
    laguerre5(f,r = 1, a = 0)
    
Performs a 10-point translated gauss-laguerre integration of the form:

```
∫exp(-ry)f(y) dy,    y ∈ (a,∞)
```
"""
function laguerre10(Base.@specialize(f),r = 1.,a = 0.)
    return _laguerrex(f,r,a,laguerre10_u,laguerre10_w)
 end
function _laguerrex(Base.@specialize(f),r,a,u,w)
    k = exp(-r*a)/r
    u1,w1 = u[1],w[1]
    rinv = 1/r
    x1 = u1*rinv + a
    f1 = f(x1)
    res = w1*f1
    l = length(u)
    for i in 2:l
        xi = u[i]*rinv + a
        res += w[i]*f(xi)
    end
    return res*k
end

"""
    evalpolyint(x0,p)

Evaluates the indefinite integral of ∑x^(k-1)*p[k] at a point x0.

"""
function evalpolyint(x,p)
    f(x) = x[2]/x[1]
    p_integral = map(f,enumerate(p))
    return evalpoly(x,p_integral)*x
end


