### functions to use:
mix_mean(pᵢ,pⱼ,k=0) = 0.5*(pᵢ+pⱼ)*(1-k)
mix_geomean(pᵢ,pⱼ,k=0) = sqrt(pᵢ*pⱼ)*(1-k)
mix_powmean(pᵢ,pⱼ,k=0,n=2) = (1-k)*(0.5*(pᵢ^n + pⱼ^n))^(1/n)
mix_mean3(pᵢ,pⱼ,k=0) = (1-k)*(0.5*(cbrt(pᵢ) + cbrt(pⱼ)))^3

##special lambda with custom k
function mix_lambda(λᵢ,λⱼ,k)
    return k + sqrt((λᵢ - k) * (λⱼ - k))
end

struct MixLambda{K}
    k::K
end

(m::MixLambda{K})(λᵢ,λⱼ,z) where K = mix_lambda(λᵢ,λⱼ,m.k)

# for use in pair_mix
mix_HudsenMcCoubrey(ϵᵢ,ϵⱼ,σᵢ,σⱼ,σᵢⱼ) = sqrt(ϵᵢ*ϵⱼ)*(σᵢ^3 * σⱼ^3)/σᵢⱼ^6
mix_HudsenMcCoubreysqrt(ϵᵢ,ϵⱼ,σᵢ,σⱼ,σᵢⱼ) = sqrt(ϵᵢ*ϵⱼ*(σᵢ^3 * σⱼ^3))/σᵢⱼ^3
mix_lambda_squarewell(λᵢ,λⱼ,σᵢ,σⱼ,σᵢⱼ) = (σᵢ*λᵢ + σⱼ*λⱼ)/(σᵢ + σⱼ)

#for use in ijab_mix
mix_ijab_elliott(βᵢ,βⱼ,σᵢ,σⱼ,σᵢⱼ) = sqrt(βᵢ * βⱼ * σᵢ^3 * σⱼ^3) / σᵢⱼ^3
mix_ijab_elliott(βᵢ,βⱼ,σᵢ,σⱼ,σᵢⱼ) = sqrt(βᵢ * βⱼ)

#throw error if the pair_mix function requires qij, but just a vector is provided.
__requires_qij(x) = false
__requires_qij(::typeof(mix_HudsenMcCoubrey)) = true
__requires_qij(::typeof(mix_HudsenMcCoubreysqrt)) = true

@noinline function __qij_error(f)
    throw(ArgumentError(lazy"$f requires a matrix of elements. A vector was provided."))
end

