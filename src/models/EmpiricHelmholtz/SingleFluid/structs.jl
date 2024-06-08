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

struct GaoBTerm
    active::Bool
    n::Vector{Float64}
    t::Vector{Float64}
    d::Vector{Float64}
    eta::Vector{Float64}
    beta::Vector{Float64}
    gamma::Vector{Float64}
    epsilon::Vector{Float64}
    b::Vector{Float64}
    function GaoBTerm(n,t,d,eta,beta,gamma,epsilon,b)
        @assert length(eta) == length(beta) == length(gamma) == length(epsilon) == length(b)
        @assert length(eta) == length(n) == length(t) == length(d)
        active = (length(n) != 0)
        return new(active,n,t,d,eta,beta,gamma,epsilon,b)
    end
end

GaoBTerm() = GaoBTerm(Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[])

struct NonAnalyticTerm 
    active::Bool
    A::Vector{Float64}
    B::Vector{Float64}
    C::Vector{Float64}
    D::Vector{Float64}
    a::Vector{Float64}
    b::Vector{Float64}
    beta::Vector{Float64}
    n::Vector{Float64}
    function NonAnalyticTerm(A,B,C,D,a,b,beta,n)
        @assert length(A) == length(B) == length(C) == length(D)
        @assert length(A) == length(a) == length(b) == length(beta)
        @assert length(beta) == length(n)
        active = (length(n) != 0)
        return new(active,A,B,C,D,a,b,beta,n)
    end
end

NonAnalyticTerm() = NonAnalyticTerm(Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[])

mutable struct Associating2BTerm #mutable because someone would want to fit this?
    active::Bool
    epsilonbar::Float64
    kappabar::Float64
    a::Float64
    m::Float64
    vbarn::Float64
    function Associating2BTerm(epsilonbar,kappabar,a,m,vbarn)
        active = (kappabar != 0.0)
        return new(active,epsilonbar,kappabar,a,m,vbarn)
    end
end
Associating2BTerm() = Associating2BTerm(0.0,0.0,0.0,0.0,0.0)

#we store power, exponential and gaussian terms inline, because those are the most used.

struct SingleFluidResidualParam <: MultiParameterParam
    iterators::Vector{UnitRange{Int}}
    n::Vector{Float64}
    t::Vector{Float64}
    d::Vector{Float64}
    l::Vector{Float64}
    g::Vector{Float64}
    eta::Vector{Float64}
    beta::Vector{Float64}
    gamma::Vector{Float64}
    epsilon::Vector{Float64}
    gao_b::GaoBTerm
    na::NonAnalyticTerm
    assoc::Associating2BTerm
    
    function SingleFluidResidualParam(n,t,d,l = Float64[],g = ones(length(l)),
        eta = Float64[],beta = Float64[],gamma = Float64[], epsilon = Float64[],
        ;gao_b = GaoBTerm(),
        na = NonAnalyticTerm(),
        assoc = Associating2BTerm())

        param = new(Vector{UnitRange{Int}}(undef,0),n,t,d,l,g,eta,beta,gamma,epsilon,gao_b,na,assoc)
        _calc_iterators!(param)
        return param
    end
end

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
        def_terms = "Lead terms: $(param.a1) + $(param.a2)*τ + $(param.c0)*log(τ)"
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

    if hasfield(typeof(param),:iterators)
        k_pol,k_exp,k_gauss = param.iterators
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

    if !__has_extra_params(param)
        show_pairs(io,res,quote_string = false)
        return nothing
    end
    #special terms
    if hasfield(typeof(param),:na)
        if param.na.active
            push!(res,"Non Analytic terms: $(length(param.na.beta))")
        end
    end

    if hasfield(typeof(param),:gao_b)
        if param.gao_b.active
            push!(res,"Gao-b terms: $(length(param.gao_b.b))")
        end
    end
    if hasfield(typeof(param),:assoc)
        if param.assoc.active
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
    resize!(param.iterators,3)
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

    function SingleFluidProperties(Mw,Tr,rhor,lb_volume,
        Tc = NaN,Pc = NaN,rhoc = NaN,
        Ttp = NaN,ptp = NaN, rhov_tp = NaN,rhol_tp = NaN,
        acentric_factor = NaN, Rgas = R̄)
        return new(Mw,Tr,rhor,lb_volume,Tc,Pc,rhoc,Ttp,ptp,rhov_tp,rhol_tp,acentric_factor,Rgas)
    end
end

is_splittable(::SingleFluidProperties) = false

const ESFProperties = SingleFluidProperties
const ESFIdealParam = SingleFluidIdealParam
const ESFResidualParam = SingleFluidResidualParam
