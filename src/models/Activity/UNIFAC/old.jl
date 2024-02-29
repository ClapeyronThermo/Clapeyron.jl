#=
function activity_coefficient(model::UNIFACModel,V,T,z)
    return exp.(@f(lnγ_comb)+ @f(lnγ_res))
end

function lnγ_comb(model::UNIFACModel,V,T,z)
    x = z ./ sum(z)
    r =model.unifac_cache.r
    q =model.unifac_cache.q
    q_p = model.unifac_cache.q_p
    Φ = r/dot(x,r)
    Φ_p = q_p/dot(x,q_p)
    θ = q/dot(x,q)
    lnγ_comb = @. log(Φ_p)+(1-Φ_p)-5*q*(log(Φ/θ)+(1-Φ/θ))
    return lnγ_comb
end

function lnγ_SG(model::UNIFACModel,V,T,z)

    x = z ./ sum(z)

    r =model.unifac_cache.r
    q =model.unifac_cache.q

    Φ = r/dot(x,r)
    θ = q/dot(x,q)
    lnγ_SG = @. -5*q*(log(Φ/θ)+(1-Φ/θ))
    return lnγ_SG
end

function lnγ_res(model::UNIFACModel,V,T,z)
    v  = model.groups.n_flattenedgroups
    _ψ = @f(Ψ)
    lnΓ_ = @f(lnΓ,_ψ)
    lnΓi_ = @f(lnΓi,_ψ)
    lnγ_res_ = [sum(v[i][k].*(lnΓ_[k].-lnΓi_[i][k]) for k ∈ @groups) for i ∈ @comps]
    return lnγ_res_
end

function lnΓ(model::UNIFACModel,V,T,z,ψ = @f(ψ))
    Q = model.params.Q.values
    v  = model.groups.n_flattenedgroups
    x = z ./ sum(z)
    X = sum(v[i][:]*x[i] for i ∈ @comps) ./ sum(sum(v[i][k]*x[i] for k ∈ @groups) for i ∈ @comps)
    θ = X.*Q / dot(X,Q)
    lnΓ_ = Q.*(1 .-log.(sum(θ[m]*ψ[m,:] for m ∈ @groups)) .- sum(θ[m]*ψ[:,m]./sum(θ[n]*ψ[n,m] for n ∈ @groups) for m ∈ @groups))
    return lnΓ_
end

function lnΓi(model::UNIFACModel,V,T,z,ψ = @f(ψ))
    Q = model.params.Q.values
    v  = model.groups.n_flattenedgroups
    ψ = @f(Ψ)
    X = [v[i][:] ./ sum(v[i][k] for k ∈ @groups) for i ∈ @comps]
    θ = [X[i][:].*Q ./ sum(X[i][n]*Q[n] for n ∈ @groups) for i ∈ @comps]
    lnΓi_ = [Q.*(1 .-log.(sum(θ[i][m]*ψ[m,:] for m ∈ @groups)) .- sum(θ[i][m]*ψ[:,m]./sum(θ[i][n]*ψ[n,m] for n ∈ @groups) for m ∈ @groups)) for i ∈ @comps]
    return lnΓi_
end

function Ψ(model::UNIFACModel,V,T,z)
    A = model.params.A.values
    B = model.params.B.values
    C = model.params.C.values
    return @. exp(-(A+B*T+C*T^2)/T)
end


=#