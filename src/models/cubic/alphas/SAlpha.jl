"""
    SAlpha{N,A} <: AlphaModel

A static wrapper over an alpha model, that returns a static vector instead of a dynamic one.

"""
struct SAlpha{N,A} <: AlphaModel
    alpha::A
end

component_list(model::SAlpha) = component_list(model.alpha)
each_split_model(model::SAlpha{N,A},I) where {N,A} = each_split_model(model.alpha,I)

function Clapeyron.α_function(model::CubicModel,V,T,z,alpha_model::SAlpha{NN, <: MathiasCopemanAlphaModel}) where {NN}
    Tc = model.params.Tc.values
    c1 = alpha_model.params.c1.values
    c2 = alpha_model.params.c2.values
    c3 = alpha_model.params.c3.values
    _1 = oneunit(eltype(alpha_model))
    _0 = Base.promote_eltype(model,T)
    n_dynamic = length(model)
    α0 = ntuple(Val(NN)) do i
        if i <= n_dynamic
            Tr = T/Tc[i]
            _1_Tr = 1 - sqrt(Tr)
            evalpoly(_1_Tr,(_1,c1[i],c2[i],c3[i]))^2
        else
            _0
        end
    end
    αvec = SVector(α0)
    α = @view(αvec,Base.OneTo(n_dynamic))
    return α
end

function Clapeyron.α_function(model::CubicModel,V,T,z,alpha_model::SAlpha{NN, <: SoaveAlphaModel}) where {NN}
    Tc = model.params.Tc.values
    ω  = alpha_model.params.acentricfactor.values
    coeff = α_m(model,alpha_model)
    _0 = Base.promote_eltype(model,T)
    n_dynamic = length(model)
    α0 = ntuple(Val(NN)) do i
        if i <= n_dynamic
            ωi = ω[i]
            Tr = T/Tc[i]
            m = evalpoly(ωi,coeff)
            (1+m*(1-sqrt(Tr)))^2
        else
            _0
        end
    end
    αvec = SVector(α0)
    α = @view(αvec,Base.OneTo(n_dynamic))
    return α
end

function Clapeyron.α_function(model::CubicModel,V,T,z,alpha_model::SAlpha{NN, <: GeneralizedSuaveAlphaModel}) where {NN}
    Tc = model.params.Tc.values
    _0 = Base.promote_eltype(model,T)
    n_dynamic = length(model)
    α0 = ntuple(Val(NN)) do i
        if i <= n_dynamic
            Tr = T/Tc[i]
            m = α_m(model,alpha_model,i)
            (1+m*(1-sqrt(Tr)))^2
        else
            _0
        end
    end
    αvec = SVector(α0)
    α = @view(αvec,Base.OneTo(n_dynamic))
    return α
end

function Clapeyron.α_function(model::CubicModel,V,T,z,alpha_model::SAlpha{NN, <: TwuAlphaModel}) where {NN}
    Tc = model.params.Tc.values
    _M  = alpha_model.params.M.values
    _N  = alpha_model.params.N.values
    _L  = alpha_model.params.L.values
    _0 = Base.promote_eltype(model,T)
    n_dynamic = length(model)
    α0 = ntuple(Val(NN)) do i
        if i <= n_dynamic
            M = _M[i]
            N = _N[i]
            L = _L[i]
            Tr = T/Tc[i]
            Tr^(N*(M-1))*exp(L*(1-Tr^(N*M)))
        else
            _0
        end
    end
    αvec = SVector(α0)
    α = @view(αvec,Base.OneTo(n_dynamic))
    return α
end

function Clapeyron.α_function(model::CubicModel,V,T,z,alpha_model::SAlpha{NN, <: RKPRAlphaModel}) where {NN}
    k1 = alpha_model.params.k1.values
    k2 = alpha_model.params.k2.values
    Tc = model.params.Tc.values
    _0 = Base.promote_eltype(model,T)
    n_dynamic = length(model)
    α0 = ntuple(Val(NN)) do i
        if i <= n_dynamic
            Tr = T/Tc[i]
            k2i = k2[i]
            k1i = k1[i]
            ((k2i + 1)/(k2i + Tr))^k1i
        else
            _0
        end
    end
    αvec = SVector(α0)
    α = @view(αvec,Base.OneTo(n_dynamic))
    return α
end

statify_alpha(αmodel,n) = αmodel
statify_alpha(αmodel::TwuAlphaModel,n) = SAlpha{n,typeof(αmodel)}(αmodel)
statify_alpha(αmodel::GeneralizedSuaveAlphaModel,n) = SAlpha{n,typeof(αmodel)}(αmodel)
statify_alpha(αmodel::SoaveAlphaModel,n) = SAlpha{n,typeof(αmodel)}(αmodel)
statify_alpha(αmodel::MathiasCopemanAlphaModel,n) = SAlpha{n,typeof(αmodel)}(αmodel)
statify_alpha(αmodel::RKPRAlphaModel,n) = SAlpha{n,typeof(αmodel)}(αmodel)
