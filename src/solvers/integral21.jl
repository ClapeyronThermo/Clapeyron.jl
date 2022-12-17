const gk21_x = [-0.9956571630258080807355,-0.973906528517171720078,-0.9301574913557082260012,-0.8650633666889845107321,-0.7808177265864168970637,-0.6794095682990244062343,-0.562757134668604683339,
        -0.4333953941292471907993,-0.294392862701460198131,-0.1488743389816312108848,0,0.1488743389816312108848,0.2943928627014601981311,0.4333953941292471907993,
        0.562757134668604683339,0.6794095682990244062343,0.7808177265864168970637,0.865063366688984510732,0.9301574913557082260012,0.973906528517171720078,0.9956571630258080807355]

const gk21_w = [0.0116946388673718742781,0.0325581623079647274788,0.0547558965743519960314,0.075039674810919952767,0.093125454583697605535,0.1093871588022976418992,0.123491976262065851078,
        0.134709217311473325928,0.142775938577060080797,0.1477391049013384913748,0.149445554002916905665,0.1477391049013384913748,0.1427759385770600807971,0.134709217311473325928,
        0.123491976262065851078,0.109387158802297641899,0.093125454583697605535,0.075039674810919952767,0.05475589657435199603138,0.032558162307964727479,0.0116946388673718742781]


"""
    integral21(f,a,b)

Performs integration via a 21-point Gauss-Kronrod rule 
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

function integral21(@specialize(f),a::ForwardDiff.Dual{T},b) where {T}
    ā,b̄ = ForwardDiff.value(a),b
    f̄ = integral21(f,ā,b̄)
    df = -f(ā)
    return ForwardDiff.dual_definition_retval(Val{T}(), f̄, df, ForwardDiff.partials(a))
end

function integral21(@specialize(f),a,b::ForwardDiff.Dual{T}) where {T}
    ā,b̄ = a,ForwardDiff.value(b)
    f̄ = integral21(f,ā,b̄)
    df = f(b̄)
    return ForwardDiff.dual_definition_retval(Val{T}(), f̄, df, ForwardDiff.partials(b))
end

function integral21(@specialize(f),a::ForwardDiff.Dual{T},b::ForwardDiff.Dual{T}) where {T}
    ā,b̄ = ForwardDiff.value(a),ForwardDiff.value(b)
    f̄ = integral21(f,ā,b̄)
    dfa = - f(ā)
    dfb = f(b̄) 
    return ForwardDiff.dual_definition_retval(Val{T}(), f̄, dfa, ForwardDiff.partials(a), dfb, ForwardDiff.partials(b))
end
