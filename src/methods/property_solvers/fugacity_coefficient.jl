# Function to compute fugacity coefficient
function lnϕ(model::EoSModel, p, T, z=SA[1.],cache = nothing;
            phase=:unknown,
            vol0=nothing,
            threaded = true,
            vol = volume(model,p,T,z;phase,vol0,threaded))


    RT = Rgas(model)*T
    logZ = log(p*vol/RT/sum(z))
    nc = length(z)
    if cache !== nothing
        result,aux,lnϕ,∂lnϕ∂n,∂lnϕ∂P,∂P∂n,∂lnϕ∂T,hconfig = cache
        if nc == 1
            lnϕ[1] = VT_lnϕ_pure(model,vol/sum(z),T,p)
        else
            aux .= 0
            aux[1] = vol
            aux[2:nc+1] = z
            gconfig = Solvers._GradientConfig(hconfig)
            F_res(model, V, T, z) = eos_res(model, V, T, z)
            fun(aux) = F_res(model, aux[1], T, @view(aux[2:nc+1]))
            _result = ForwardDiff.gradient!(result, fun, aux, gconfig, Val{false}())
            dresult = DiffResults.gradient(_result)
            μ_res = @view dresult[2:nc+1]
            lnϕ .= μ_res ./ RT .- logZ
        end
    else
        μ_res = VT_chemical_potential_res(model, vol, T, z)
        lnϕ = μ_res/RT .- logZ
    end
    return lnϕ, vol
end

function lnϕ(model::IdealModel, p, T, z=SA[1.],cache = nothing;
            phase=:unknown,
            vol0=nothing,
            threaded = true,
            vol = volume(model,p,T,z;phase,vol0,threaded))

    lnϕ = FillArrays.Zeros(length(z))
    return lnϕ, vol
end

function lnϕ!(cache::Tuple, model::EoSModel, p, T, z=SA[1.];
            phase=:unknown,
            vol0=nothing,
            threaded = true,
            vol = volume(model,p,T,z;phase,vol0,threaded))

    model isa IdealModel && return (fill!(cache[3],0.0),vol)

    return lnϕ(model,p,T,z,cache;vol)
end

function lnϕ!(lnϕ::AbstractVector, model::EoSModel, p, T, z=SA[1.],cache = nothing;
            phase=:unknown,
            vol0=nothing,
            threaded = true,
            vol = volume(model,p,T,z;phase,vol0,threaded))

    if model isa IdealModel
        lnϕ .= 0
        return lnϕ,vol
    end

    if cache != nothing
        lnϕ(model,p,T,z,cache;vol)
        result,aux,lnϕw,∂lnϕ∂n,∂lnϕ∂P,∂P∂n,∂lnϕ∂T,hconfig = cache
        vol = aux[1]
        lnϕ .= lnϕw
    else
        RT = Rgas(model)*T
        μ_res = VT_chemical_potential_res!(lnϕ,model, vol, T, z)
        Z = p*vol/RT/sum(z)
        lnϕ .= μ_res ./ RT .- log(Z)
    end
    return lnϕ, vol
end

function VT_lnϕ_pure(model,V,T,p = pressure(model,V,T))
    RT = Rgas(model)*T
    p_res = p - RT/V
    μ_res = eos_res(model,V,T) + p_res*V
    Z = p*V/RT
    return μ_res/RT - log(Z)
end

function ∑zlogϕ(model::EoSModel, p, T, z=SA[1.],cache = nothing;
            phase=:unknown,
            vol0=nothing,
            threaded = true,
            vol = volume(model,p,T,z;phase,vol0,threaded))

    return VT_∑zlogϕ(model,vol,T,z),vol
end

function VT_∑zlogϕ(model,V,T,z)
    RT = Rgas(model)*T
    n = sum(z)
    A, ∂A∂V, ∂A∂T = ∂f_res_vec(model,V,T,z)
    PrV = ifelse(iszero(1/V),zero(∂A∂V),- V*∂A∂V)
    g_res = A + PrV
    logZ = log1p(∂A∂V*V/(n*RT))
    ∑zlogϕi = g_res/RT - n*logZ
    return ∑zlogϕi
end

struct ∂lnϕTag end

function ∂lnϕ_cache(model::EoSModel, p, T, z, ::Val{B}) where B
    TT = Base.promote_eltype(model,p,T,z)
    lnϕ = zeros(TT,length(model))
    aux = zeros(TT,length(model) + 1 + B)
    ∂lnϕ∂n = lnϕ * transpose(lnϕ)
    result = DiffResults.HessianResult(aux)
    ∂lnϕ∂n = lnϕ * transpose(lnϕ)
    ∂lnϕ∂P = similar(lnϕ)
    ∂P∂n = similar(lnϕ)
    hconfig = ForwardDiff.HessianConfig((∂lnϕTag(),model,p,T,z),result,aux)
    ∂tag = ForwardDiff.tagtype(hconfig.jacobian_config.duals[1][1])
    if has_lnγ_impl(__γ_unwrap(model))
        jcache = similar(aux) 
        dlnγdT_cache = Vector{ForwardDiff.Dual{∂tag,TT,1}}(undef,length(model))
    else
        jcache = aux
        dlnγdT_cache = Vector{ForwardDiff.Dual{∂tag,TT,1}}(undef,0)
    end
    if B
        ∂lnϕ∂T = similar(lnϕ)
    else
        ∂lnϕ∂T = lnϕ
    end
    result,aux,lnϕ,∂lnϕ∂n,∂lnϕ∂P,∂P∂n,∂lnϕ∂T,hconfig,jcache,dlnγdT_cache
end

# Function to compute fugacity coefficient and its pressure and composition derivatives
function ∂lnϕ∂n∂P(model::EoSModel, p, T, z=SA[1.], cache = ∂lnϕ_cache(model,p,T,z,Val{false}());
            phase=:unknown,
            vol0=nothing,
            threaded = true,
            vol = volume(model,p,T,z;phase,vol0,threaded))

    result,aux,lnϕ,∂lnϕ∂n,∂lnϕ∂P,∂P∂n,∂lnϕ∂T,hconfig = cache
    RT = Rgas(model)*T
    V = vol
    n = sum(z)
    Z = p*V/RT/n

    F_res(model, V, T, z) = eos_res(model, V, T, z) / RT
    fun(aux) = F_res(model, aux[1], T, @view(aux[2:(ncomponents+1)]))

    ncomponents = length(z)
    aux[1] = V
    aux[2:end] = z
    result = ForwardDiff.hessian!(result, fun, aux, hconfig, Val{false}())

    F = DiffResults.value(result)
    ∂F = DiffResults.gradient(result)
    ∂2F = DiffResults.hessian(result)

    ∂F∂V = ∂F[1]
    ∂F∂n = @view ∂F[2:(ncomponents+1)]

    ∂2F∂V2 = ∂2F[1, 1]
    ∂2F∂n2 = @view ∂2F[2:(ncomponents+1), 2:(ncomponents+1)]
    ∂2F∂n∂V = @view ∂2F[1, 2:(ncomponents+1)]

    ∂P∂V = -RT*∂2F∂V2 - n*RT/V^2
    ∂P∂n .= - RT .* ∂2F∂n∂V .+ RT ./ V
    lnϕ .= ∂F∂n .- log(Z)
    for i in 1:ncomponents
        ∂P∂ni = ∂P∂n[i]
        ∂V∂ni = - ∂P∂ni/∂P∂V
        ∂lnϕ∂P[i] = ∂V∂ni/RT - 1/p
        for j in 1:ncomponents
            ∂lnϕ∂n[i,j] = ∂2F∂n2[i,j] + 1/n + (∂P∂ni * ∂P∂n[j])/∂P∂V/RT
        end
    end
    return lnϕ, ∂lnϕ∂n, ∂lnϕ∂P, V
end

# Function to compute fugacity coefficient pressure derivatives
function ∂lnϕ∂P(model::EoSModel, p, T, z=SA[1.], cache = ∂lnϕ_cache(model,p,T,z,Val{false}());
            phase=:unknown,
            vol0=nothing,
            threaded = true,
            vol = volume(model,p,T,z;phase,vol0,threaded))

    RT = Rgas(model)*T
    nc = length(z)
    if cache === nothing
        f1(w) = pressure_res(model,vol,T,w)
        ∇pᵣ = ForwardDiff.gradient(f1,z)
        ∂lnϕ∂P = ∇pᵣ
        _,∂p∂V = p∂p∂V(model,vol,T,z)
        ∂V∂p = -1/∂p∂V
        ∂lnϕ∂P .= ∂V∂p .* ∇pᵣ ./ RT .- (1/p) .+ (∂V∂p/vol)
        return ∂lnϕ∂P,vol
    end
    
    result,aux,lnϕ,∂lnϕ∂n,∂lnϕ∂P,∂P∂n,∂lnϕ∂T,hconfig = cache
    aux .= 0
    aux[nc+1] = vol
    aux[1:nc] = z
    gconfig = Solvers._GradientConfig(hconfig)
    fun(aux) = pressure_res(model,aux[nc+1],T,aux[1:nc])
    _result = ForwardDiff.gradient!(result, fun, aux, gconfig, Val{false}())
    dresult = DiffResults.gradient(_result)
    ∇pᵣ = @view dresult[1:nc]
    ∂pᵣ∂V = dresult[nc+1]
    ∂p∂V = ∂pᵣ∂V - sum(z)*RT/(vol*vol)
    ∂V∂p = -1/∂p∂V
    ∂lnϕ∂P .= ∂V∂p .* ∇pᵣ ./ RT .- (1/p) .+ (∂V∂p/vol)
    return ∂lnϕ∂P,vol
end

function dardT(model,v,T,z)
    f(_T) = eos_res(model,v,_T,z)/(Rgas(model)*_T)
    return Solvers.derivative(f,T)
end

# Function to compute fugacity coefficient temperature derivatives
function ∂lnϕ∂T(model::EoSModel, p, T, z=SA[1.], cache = ∂lnϕ_cache(model,p,T,z,Val{false}());
            phase=:unknown,
            vol0=nothing,
            threaded = true,
            vol = volume(model,p,T,z;phase,vol0,threaded))

    RT = Rgas(model)*T
    fun_T(w) = dardT(model,vol,T,w)
    fun_v(w) = -pressure_res(model,vol,T,w)/RT
    if cache === nothing
        ∇sᵣ = ForwardDiff.gradient(fun_T,z)
        ∇pᵣ = ForwardDiff.gradient(fun_v,z)
        d2A = f_hess(model,vol,T,z)
        ∂p∂V = -d2A[1,1]
        ∂p∂T = -d2A[1,2]
        ∂V∂T = -∂p∂T/∂p∂V
        ∂lnZ∂T = ∂V∂T/vol - 1/T
        return ∇sᵣ .- ∇pᵣ .* ∂p∂T ./ ∂p∂V .- ∂lnZ∂T,vol
    else
        n = sum(z)
        nc = length(z)
        result,aux,lnϕ,∂lnϕ∂n,∇pᵣ,∇sᵣ,∂lnϕ∂T,hconfig = cache
        aux .= 0
        aux[1:nc] = z
        aux[nc+1] = vol
        aux[nc+2] = T
        fun_v2(aux) = -pressure_res(model,aux[nc+1],aux[nc+2],@view(aux[1:nc]))/(Rgas(model)*aux[nc+2])
        fun_t2(aux) = dardT(model,vol,T,@view(aux[1:nc]))
        gconfig = Solvers._GradientConfig(hconfig)
        _result_p = ForwardDiff.gradient!(result, fun_v2, aux, gconfig, Val{false}())
        dresult_p = DiffResults.gradient(_result_p)
        ∇pᵣ .= @view dresult_p[1:nc]
        ∂pᵣ∂V = dresult_p[nc+1]
        #@show ∂pᵣ∂V
        ∂pᵣ∂T = dresult_p[nc+2]
        ∂p∂V = RT*(-∂pᵣ∂V - n/(vol*vol))
        ∂p∂T = -RT*∂pᵣ∂T + p/T
        ∂V∂T = -∂p∂T/∂p∂V
        ∂lnZ∂T = ∂V∂T/vol - 1/T
        _result_s = ForwardDiff.gradient!(result, fun_t2, aux, gconfig, Val{false}())
        dresult_s = DiffResults.gradient(_result_s)
        ∇sᵣ .= @view dresult_s[1:nc]
        ∂lnϕ∂T .= ∇sᵣ .- ∇pᵣ .* ∂p∂T ./ ∂p∂V .- ∂lnZ∂T
        return ∂lnϕ∂T,vol
    end
end

# Function to compute fugacity coefficient and its temperature, pressure and composition derivatives
function ∂lnϕ∂n∂P∂T(model::EoSModel, p, T, z=SA[1.],cache = ∂lnϕ_cache(model,p,T,z,Val{true}());
            phase=:unknown,
            vol0=nothing,
            threaded = true,
            vol = volume(model,p,T,z;phase,vol0,threaded))

    result,aux,lnϕ,∂lnϕ∂n,∂lnϕ∂P,∂P∂n,∂lnϕ∂T,hconfig = cache
    RT = Rgas(model)*T
    V = vol
    n = sum(z)
    Z = p*V/RT/n

    ncomponents = length(z)
    aux[1] = V
    aux[2] = T
    aux[3:end] .= z
    F_res(model, V, T, z) = eos_res(model, V, T, z) / (Rgas(model)*T)
    fun(aux) = F_res(model, aux[1], aux[2], @view(aux[3:(ncomponents+2)]))
    result = ForwardDiff.hessian!(result, fun, aux, hconfig, Val{false}())

    F = DiffResults.value(result)
    ∂F = DiffResults.gradient(result)
    ∂2F = DiffResults.hessian(result)
    ∂F∂V = ∂F[1]
    ∂F∂T = ∂F[2]
    ∂F∂n = @view ∂F[3:(ncomponents+2)]

    ∂2F∂V2 = ∂2F[1, 1]
    ∂2F∂n2 = @view ∂2F[3:(ncomponents+2), 3:(ncomponents+2)]
    ∂2F∂T2 = ∂2F[2, 2]
    ∂2F∂n∂V = @view ∂2F[1, 3:(ncomponents+2)]
    ∂2F∂n∂T = @view ∂2F[2, 3:(ncomponents+2)]
    ∂2F∂V∂T = ∂2F[1, 2]
    lnϕ .= ∂F∂n .- log(Z)
    ∂P∂V = -RT*∂2F∂V2 - n*RT/V^2
    ∂P∂n .= -RT .* ∂2F∂n∂V .+ RT ./ V
    ∂P∂T = -RT*∂2F∂V∂T + p/T
    for i in 1:ncomponents
        ∂P∂ni = ∂P∂n[i]
        ∂V∂ni = - ∂P∂ni/∂P∂V
        ∂lnϕ∂P[i] = ∂V∂ni/RT - 1/p
        ∂lnϕ∂T[i] = ∂2F∂n∂T[i] + 1/T - (∂V∂ni*∂P∂T)/RT
        for j in 1:ncomponents
            ∂lnϕ∂n[i,j] = ∂2F∂n2[i,j] + 1/n + (∂P∂ni * ∂P∂n[j])/∂P∂V/RT
        end
    end

    return lnϕ, ∂lnϕ∂n, ∂lnϕ∂P, ∂lnϕ∂T, V
end

#functions used to overload gamma-phi models
function modified_lnϕ(model, p, T, z, cache; phase = :unknown, vol0 = nothing)
    lnϕz,vz = lnϕ(model, p, T, z, cache; phase, vol0)
    if isnan(vz)
        lnϕz,vz = lnϕ(model, p, T, z, cache; phase)
    end
    return lnϕz,vz
end

function modified_∂lnϕ∂n(model, p, T, z, cache; phase = :unknown, vol0 = nothing)
    lnϕ, ∂lnϕ∂n, _, vol = ∂lnϕ∂n∂P(model, p, T, z, cache; phase, vol0)
    return lnϕ,∂lnϕ∂n,vol
end

#=
VT-based versions

instead of working in terms of lnphi, we work directly in terms of lnf = ln(phi*p)

 Z = p*V/RT/sum(z)
lnϕ .= μ_res ./ RT .- log(Z)
ln(f) = ln(ϕ*p) = log(ϕ) + log(p)
ln(f) = μ_res ./ RT .- log(V/RT/sum(z))
- logZ - log(p) =
=#

function lnf(model, V, T, z,cache = nothing)
    RT = Rgas(model)*T
    n = sum(z)
    logZp = log(V/RT/n)
    nc = length(z)
    TT = Base.promote_eltype(model,V,T,z)
    F_res(_model, _V, _T, _z) = eos_res(_model, _V, _T, _z)
    fun(aux) = F_res(model, aux[1], T, @view(aux[2:nc+1]))
    if cache !== nothing
        result,aux,lnf,∂lnϕ∂n,∂lnϕ∂P,∂P∂n,∂lnϕ∂T,hconfig = cache
        if nc == 1
            lnf1,p = VT_lnf_pure(model,V,T)
            lnf[1] = lnf1
        else
            aux .= 0
            aux[1] = V
            aux[2:nc+1] = z
            gconfig = Solvers._GradientConfig(hconfig)
            _result = ForwardDiff.gradient!(result, fun, aux, gconfig, Val{false}())
            dresult = DiffResults.gradient(_result)
            dfdv = dresult[1]
            μ_res = @view dresult[2:nc+1]
            p = -(dfdv - RT*n/V)
            lnf .= μ_res ./ RT .- logZp
        end
    else
        μ_res = VT_chemical_potential_res(model,V,T,z)
        p = pressure(model,V,T,z)
        lnf = μ_res ./ RT .- logZp
    end
    return lnf, p
end

function VT_lnf_pure(model,V,T)
    RT = Rgas(model)*T
    f(dV) = eos_res(model,dV,T,SA[1.0])
    F,dFdV = Solvers.f∂f(f,V)
    p_res = -dFdV
    μ_res = eos_res(model,V,T,SA[1.0]) + p_res*V
    Zp = V/RT
    lnf = μ_res/RT - log(Zp)
    p = p_res + RT/V
    return lnf,p
end

function ∂lnf∂n∂V(model, V, T, z, cache = ∂lnϕ_cache(model,V,T,z,Val{false}()))
    
    result,aux,lnf,∂lnf∂n,∂lnf∂V,∂P∂n,∂lnf∂T,hconfig = cache
    
    RT = Rgas(model)*T
    n = sum(z)
    Zp = V/RT/n
    ncomponents = length(z)
    
    F_res(model, V, T, z) = eos_res(model, V, T, z) / RT
    fun(aux) = F_res(model, aux[1], T, @view(aux[2:(ncomponents+1)]))
  
    aux[1] = V
    aux[2:end] = z
    result = ForwardDiff.hessian!(result, fun, aux, hconfig, Val{false}())
    
    F = DiffResults.value(result)
    ∂F = DiffResults.gradient(result)
    ∂2F = DiffResults.hessian(result)
    
    ∂F∂V = ∂F[1]
    ∂F∂n = @view ∂F[2:(ncomponents+1)]

    ∂2F∂V2 = ∂2F[1, 1]
    ∂2F∂n2 = @view ∂2F[2:(ncomponents+1), 2:(ncomponents+1)]
    ∂2F∂n∂V = @view ∂2F[1, 2:(ncomponents+1)]
    p = -RT*(∂F∂V - n/V)
    ∂P∂V = -RT*∂2F∂V2 - n*RT/V^2
    ∂P∂n .= - RT .* ∂2F∂n∂V .+ RT ./ V
    lnf .= ∂F∂n .- log(Zp)
    
    for i in 1:ncomponents
        ∂P∂ni = ∂P∂n[i]
        ∂lnf∂V[i] = ∂2F∂n∂V[i] - 1/V
        for j in 1:ncomponents
            ∂lnf∂n[i,j] = ∂2F∂n2[i,j] + 1/n
        end
    end
    return lnf, ∂lnf∂n, ∂lnf∂V, ∂P∂n, p, ∂P∂V
end

function ∂lnf∂n∂V∂T(model, V, T, z, cache = ∂lnϕ_cache(model,V,T,z,Val{true}()))
    
    result,aux,lnf,∂lnf∂n,∂lnf∂V,∂P∂n,∂lnf∂T,hconfig = cache
    
    RT = Rgas(model)*T
    n = sum(z)
    Zp = V/RT/n
    ncomponents = length(z)
    
    aux[1] = V
    aux[2] = T
    aux[3:end] .= z
    F_res(model, V, T, z) = eos_res(model, V, T, z) / (Rgas(model)*T)
    fun(aux) = F_res(model, aux[1], aux[2], @view(aux[3:(ncomponents+2)]))
    result = ForwardDiff.hessian!(result, fun, aux, hconfig, Val{false}())

    F = DiffResults.value(result)
    ∂F = DiffResults.gradient(result)
    ∂2F = DiffResults.hessian(result)
    ∂F∂V = ∂F[1]
    ∂F∂T = ∂F[2]
    ∂F∂n = @view ∂F[3:(ncomponents+2)]

    ∂2F∂V2 = ∂2F[1, 1]
    ∂2F∂n2 = @view ∂2F[3:(ncomponents+2), 3:(ncomponents+2)]
    ∂2F∂T2 = ∂2F[2, 2]
    ∂2F∂n∂V = @view ∂2F[1, 3:(ncomponents+2)]
    ∂2F∂n∂T = @view ∂2F[2, 3:(ncomponents+2)]
    ∂2F∂V∂T = ∂2F[1, 2]
    
    p = -RT*(∂F∂V - n/V)
    ∂P∂V = -RT*∂2F∂V2 - n*RT/V^2
    ∂P∂n .= - RT .* ∂2F∂n∂V .+ RT ./ V
    lnf .= ∂F∂n .- log(Zp)
    ∂P∂T = -RT*∂2F∂V∂T + p/T
    for i in 1:ncomponents
        ∂P∂ni = ∂P∂n[i]
        ∂lnf∂V[i] = ∂2F∂n∂V[i] - 1/V
        ∂lnf∂T[i] = ∂2F∂n∂T[i] + 1/T
        for j in 1:ncomponents
            ∂lnf∂n[i,j] = ∂2F∂n2[i,j] + 1/n
        end
    end
    return lnf, ∂lnf∂n, ∂lnf∂V, ∂lnf∂T, ∂P∂n, p, ∂P∂V, ∂P∂T
end

function ∂lnf∂V(model::EoSModel, V, T, z=SA[1.], cache = ∂lnϕ_cache(model,V,T,z,Val{false}()))

    RT = Rgas(model)*T
    nc = length(z)
    n = sum(z)
    result,aux,_,_,∂lnf∂V,∂P∂n,_,hconfig = cache
    aux .= 0
    aux[nc+1] = V
    aux[1:nc] = z
    gconfig = Solvers._GradientConfig(hconfig)
    fun(aux) = pressure_res(model,aux[nc+1],T,aux[1:nc])
    _result = ForwardDiff.gradient!(result, fun, aux, gconfig, Val{false}())
    dresult = DiffResults.gradient(_result)
    pᵣ = DiffResults.value(_result)
    p = pᵣ + n*RT/V
    ∇pᵣ = @view dresult[1:nc]
    ∂pᵣ∂V = dresult[nc+1]
    ∂p∂V = ∂pᵣ∂V - n*RT/(V*V)
    ∂V∂p = -1/∂p∂V
    ∂lnf∂V .= - ∇pᵣ ./ RT .- 1/V
    ∂P∂n .= ∇pᵣ .+ RT ./ V
    return ∂lnf∂V, ∂P∂n, p, ∂p∂V
end