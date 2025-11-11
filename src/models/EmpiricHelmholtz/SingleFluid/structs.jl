abstract type MultiParameterParam <: EoSParam end

struct SingleFluidIdealParam <:MultiParameterParam
    a1::Float64 #a1
    a2::Float64 #a2*τ
    c0::Float64 #c0*log(τ)
    ref_a::MVector{2,Float64} #equivalent to a1/a2,but mutable and used to set reference states
    #gpe terms (Generalized Plank-Einstein)
    n_gpe::Vector{Float64} 
    t_gpe::Vector{Float64}
    c_gpe::Vector{Float64}
    d_gpe::Vector{Float64}
    
    #power terms ni*τ^ti
    n_p::Vector{Float64} 
    t_p::Vector{Float64}
    
    #GERG-2008 terms. while teorethically, they can be converted into power terms, in practice, thanks to the existence
    #of LogExpFunctions.logcosh and LogExpFunctions.logabssinh, we save evaluations.
    n_gerg::Vector{Float64}
    v_gerg::Vector{Float64}
    R0::Float64

    #c1 term: c1*τ*log(τ) . appears in IdealGasHelmholtzCP0PolyT when t = -1
    c1::Float64
    
    function SingleFluidIdealParam(a1,a2,c0,n = Float64[],t = Float64[],c = fill(1.0,length(n)),d = fill(-1.0,length(n)),n_p = Float64[], t_p  = Float64[],n_gerg = Float64[],v_gerg = Float64[],R0 = 0.0,c1 = 0.0)
        @assert length(n_gerg) == length(v_gerg)
        @assert length(n) == length(t) == length(c) == length(d)
        @assert length(n_p) == length(t_p)
        return new(a1,a2,c0,MVector(0.0,0.0),n,t,c,d,n_p,t_p,n_gerg,v_gerg,R0,c1)
    end
end

is_splittable(::SingleFluidIdealParam) = false

Base.@kwdef struct SingleFluidResidualParam <: MultiParameterParam
    polexpgauss::PolExpGaussTerm = PolExpGaussTerm()
    exp2::DoubleExpTerm = DoubleExpTerm()
    gao_b::GaoBTerm = GaoBTerm()
    na::NonAnalyticTerm = NonAnalyticTerm()
    assoc::Associating2BTerm = Associating2BTerm()
end

struct EmpiricDepartureValues <: MultiParameterParam
    polexpgauss::PolExpGaussTerm
    F::Float64
end

_calc_iterators!(m::SingleFluidResidualParam) = _calc_iterators!(m.default)
is_splittable(::SingleFluidResidualParam) = false
__has_extra_params(x::MultiParameterParam) = __has_extra_params(typeof(x))
__has_extra_params(x) = false
__has_extra_params(ℙ::Type{SingleFluidResidualParam}) = true
__has_extra_params(ℙ::Type{SingleFluidIdealParam}) = false

__type_string(x::MultiParameterParam) = __type_string(typeof(x))
__type_string(::Type{T}) where T <:MultiParameterParam = "Struct"
__type_string(ℙ::Type{SingleFluidResidualParam}) = "Residual"
__type_string(ℙ::Type{SingleFluidIdealParam}) = "Ideal"

function show_multiparameter_coeffs(io,param::MultiParameterParam)
    res = String[]

    if hasfield(typeof(param),:a1) && hasfield(typeof(param),:a2) && hasfield(typeof(param),:c0)
        def_terms1 = "Lead terms: $(param.a1) + $(param.a2)*τ"
        def_terms2 = "$(abs(param.c0))*log(τ)"
        if param.c0 > 0
            def_terms = def_terms1 * " + " * def_terms2
        else
            def_terms = def_terms1 * " - " * def_terms2
        end
        if hasfield(typeof(param),:c1)
            if param.c1 != 0
                def_terms = def_terms * " + $(param.c1)*τ*log(τ)" 
            end
        end
        push!(res,def_terms)
    end

    if hasfield(typeof(param),:n_gpe)
        if length(param.n_gpe) != 0
            push!(res,"Plank-Einstein terms: $(length(param.n_gpe))")
        end
    end

    if hasfield(typeof(param),:n_p)
        if length(param.n_p) != 0
            push!(res,"Ideal Exponential terms: $(length(param.n_p))")
        end
    end

    if hasfield(typeof(param),:n_gerg)
        if length(param.n_gerg) != 0
            push!(res,"GERG-2004 ideal terms: $(length(param.n_gerg))")
        end
    end
    if hasfield(typeof(param),:F)
        push!(res,"Fij: $(param.F)")
    end

    if hasfield(typeof(param),:polexpgauss)
        polexpgauss = param.polexpgauss
        if active_term(polexpgauss)
            k_pol,k_exp,k_gauss = polexpgauss.iterators
            if length(k_pol) != 0
                push!(res,"Polynomial power terms: $(length(k_pol))")
            end
            if length(k_exp) != 0
                push!(res,"Exponential terms: $(length(k_exp))")
            end
            if length(k_gauss) != 0
                push!(res,"Gaussian bell-shaped terms: $(length(k_gauss))")
            end
        end
    end

    #special terms

    if hasfield(typeof(param),:exp2)
        if active_term(param.exp2)
            push!(res,"Double Exponential terms: $(length(param.exp2.n))")
        end
    end

    if hasfield(typeof(param),:na)
        if active_term(param.na)
            push!(res,"Non Analytic terms: $(length(param.na.beta))")
        end
    end

    if hasfield(typeof(param),:gao_b)
        if active_term(param.gao_b)
            push!(res,"Gao-b terms: $(length(param.gao_b.b))")
        end
    end

    if hasfield(typeof(param),:assoc)
        if active_term(param.assoc)
            push!(res,"Associating terms: $(length(param.assoc.a))")
        end
    end

    show_pairs(io,res,quote_string = false)
end

function Base.show(io::IO,::MIME"text/plain",param::MultiParameterParam)
    println(io,__type_string(param)," MultiParameter coefficients:")
    show_multiparameter_coeffs(io,param)
end

function _calc_iterators!(param)
    n,t,d,l = param.n,param.t,param.d,param.l
    g = param.g
    eta,beta,gamma,epsilon = param.eta,param.beta,param.gamma,param.epsilon

    @assert length(n) == length(t) == length(d)
    @assert length(l) <= length(d)
    @assert length(g) == length(l)
    @assert length(eta) == length(beta) == length(gamma) == length(epsilon)
    #we start from the assoc term, backwards
    length_n = length(n)
    length_beta = length(beta)

    length_pol = length_n - length_beta - length(l)
    length_exp = length_n - length_beta
    length_gauss = length_n
    k_pol = 1:length_pol
    k_exp = (length_pol+1):length_exp
    k_gauss = (length_exp+1):length_gauss
    param.iterators .= (k_pol,k_exp,k_gauss)
    return param
end

struct SingleFluidProperties <: EoSParam
    #those 4 properties are the most important here:
    Mw::Float64 #Molecular Weight, g/mol
    Tr::Float64 #reducing temperature, K
    rhor::Float64 #reducing density, mol/m3
    lb_volume::Float64 #lower bound volume, mol/m3

    #the rest is optional, but recomended.
    Tc::Float64 #Critical temperature, K
    Pc::Float64 #Critical Pressure,Pa
    rhoc::Float64 #Critical density, mol/m3
    Ttp::Float64 #triple point temperature, K
    ptp::Float64 #triple point pressure, Pa
    rhov_tp::Float64 #triple point vapor volume, mol/m3
    rhol_tp::Float64 #triple point liquid volume, mol/m3
    acentricfactor::Float64 #acentric factor
    Rgas::Float64 #gas constant used
    pseudo_pure::Bool #support for pseudo-pures.

    function SingleFluidProperties(Mw,Tr,rhor,lb_volume,
        Tc = NaN,Pc = NaN,rhoc = NaN,
        Ttp = NaN,ptp = NaN, rhov_tp = NaN,rhol_tp = NaN,
        acentric_factor = NaN, Rgas = R̄,pseudo_pure = false)
        return new(Mw,Tr,rhor,lb_volume,Tc,Pc,rhoc,Ttp,ptp,rhov_tp,rhol_tp,acentric_factor,Rgas,pseudo_pure)
    end
end

is_splittable(::SingleFluidProperties) = false
Base.show(io::IO,props::SingleFluidProperties) = show_as_namedtuple(io,props)
Base.show(io::IO,::MIME"text/plain",props::SingleFluidProperties) = show_as_namedtuple(io,props)

const ESFProperties = SingleFluidProperties
const ESFIdealParam = SingleFluidIdealParam
const ESFResidualParam = SingleFluidResidualParam
