#check http://www.coolprop.org/_static/doxygen/html/class_cool_prop_1_1_saturation_ancillary_function.html

struct GenericAncEvaluator
    n::Vector{Float64}
    t::Vector{Float64}
    input_r::Float64 #reducer over input
    output_r::Float64 #reducer over output
    type::Symbol #type of evaluator
    using_input_r::Bool #on certain types,a input/input_r is used.
    superanc::Solvers.ChebyshevRangeV64 #used if EoSSuperancillaries.jl is loaded
end

function GenericAncEvaluator(n,T,input_r,output_r,type,using_input_r)
    return GenericAncEvaluator(n,T,input_r,output_r,type,using_input_r,Solvers.ChebyshevRange(Float64[],Vector{Float64}[]))
end

function _eval_generic_anc(data::GenericAncEvaluator,input)
    type = data.type
    input_r = data.input_r
    output_r = data.output_r
    xr = input/input_r
    use_xr = data.using_input_r
    xr > 1 && return zero(xr)/zero(xr)
    n = data.n
    v = data.t
    θ = (input_r-input)/input_r
    if type == :exp
        ∑nθt = evalexppoly(θ,n,v)
        if use_xr
            ∑nθt /= xr
        end
        return output_r*exp(∑nθt)
    elseif type == :noexp
        ∑nθt = evalexppoly(θ,n,v)
        return output_r*(∑nθt + 1)
    elseif type == :rational
        ∑a = evalpoly(input,n)
        ∑b = evalpoly(input,v)
        return ∑a/∑b
    elseif type == :superanc
        return Solvers.cheb_eval(data.superanc,xr)
    else
        throw(error("unrecognized type: " * string(type)))
    end
end
