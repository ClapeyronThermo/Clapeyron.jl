### functions to use:
mix_mean(p_i,p_j,k=0) = 0.5*(p_i+p_j)*(1-k)
mix_geomean(p_i,p_j,k=0) = sqrt(p_i*p_j)*(1-k)
mix_powmean(p_i,p_j,k=0,n=2) = (1-k)*(0.5*(p_i^n + p_j^n))^(1/n)
mix_mean3(p_i,p_j,k=0) = (1-k)*(0.5*(cbrt(p_i) + cbrt(p_j)))^3
##special lambda with custom k
function mix_lambda(λ_i,λ_j,k)
    return k + sqrt((λ_i - k) * (λ_j - k))
end

struct MixLambda{K}
    k::K
end

(m::MixLambda{K})(λ_i,λ_j,z) where K = mix_lambda(λ_i,λ_j,m.k)

# for use in pair_mix
mix_HudsenMcCoubrey(ϵᵢ,ϵⱼ,σᵢ,σⱼ,σᵢⱼ) = √(ϵᵢ*ϵⱼ)*(σᵢ^3 * σⱼ^3)/σᵢⱼ^6
mix_HudsenMcCoubreysqrt(ϵᵢ,ϵⱼ,σᵢ,σⱼ,σᵢⱼ) = √(ϵᵢ*ϵⱼ*(σᵢ^3 * σⱼ^3))/σᵢⱼ^3
mix_lambda_squarewell(λᵢ,λⱼ,σᵢ,σⱼ,σᵢⱼ) = (σᵢ*λᵢ + σⱼ*λⱼ)/(σᵢ + σⱼ)

#throw error if the pair_mix function requires qij, but just a vector is provided.
__requires_qij(x) = false
__requires_qij(::typeof(mix_HudsenMcCoubrey)) = true
__requires_qij(::typeof(mix_HudsenMcCoubreysqrt)) = true

@noinline function __qij_error(f)
    throw(ArgumentError(lazy"$f requires a matrix of elements. A vector was provided."))
end

