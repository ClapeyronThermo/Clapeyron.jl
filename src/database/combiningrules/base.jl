### functions to use:
mix_mean(p_i,p_j,k=0) = 0.5*(p_i+p_j)*(1-k)  
mix_geomean(p_i,p_j,k=0) = sqrt(p_i*p_j)*(1-k) 
mix_powmean(p_i,p_j,k=0,n=2) =(1-k)*(0.5*(p_i^n + p_j^n))^(1/n)

##special lambda with custom k
function mix_lambda(λ_i,λ_j,k)
    return k + sqrt((λ_i - k) * (λ_j - k))
end
# for use in pair_mix
mix_HudsenMcCoubrey(ϵᵢ,ϵⱼ,σᵢ,σⱼ,σᵢⱼ) = √(ϵᵢ*ϵⱼ)*(σᵢ^3 * σⱼ^3)/σᵢⱼ^6 
mix_lambda_squarewell(λᵢ,λⱼ,σᵢ,σⱼ,σᵢⱼ) = (σᵢ*λᵢ + σⱼ*λⱼ)/(σᵢ + σⱼ)
