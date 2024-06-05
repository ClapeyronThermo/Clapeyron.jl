# Function to accelerate Succesive Substitution by DEM method (1 eigenvalue)
function dem(xₙ, xₙ₋₁, xₙ₋₂)
    return dem!(similar(xₙ),xₙ, xₙ₋₁, xₙ₋₂)
end

function dem!(x_dem,xₙ, xₙ₋₁, xₙ₋₂,cache = (similar(x),similar(x)))
    Δxₙ₋₁,Δxₙ = cache
    Δxₙ₋₁ .= xₙ₋₁ .- xₙ₋₂
    Δxₙ .= xₙ .- xₙ₋₁
    λ = dot(Δxₙ, Δxₙ) / dot(Δxₙ, Δxₙ₋₁)
    x_dem .= xₙ .+ Δxₙ .* λ ./(1. .- λ)
end

function ass_matmul!(f::F,out,K,x) where F
    n = length(x)
    i_solved = 0
    for i in 1:n
        kxi = zero(eltype(x))
        #strategy 3
        Ki = @view K[i,:]
        @inbounds for j in 1:i_solved
            kxi += out[j]*Ki[j]
        end
        @inbounds for j in (i_solved+1):n
            kxi += x[j]*Ki[j]
        end
        out[i] = f(kxi,x[i])
        i_solved += 1
    end
end

# Function to accelerate Succesive Substitution by GDEM method (2 eigenvalues)
#=
function gdem2(xₙ, xₙ₋₁, xₙ₋₂, xₙ₋₃)
    Δxₙ₋₂ = xₙ .- xₙ₋₃
    Δxₙ₋₁ = xₙ .- xₙ₋₂
    Δxₙ = xₙ .- xₙ₋₁
    b₀₁ = dot(Δxₙ, Δxₙ₋₁)
    b₀₂ = dot(Δxₙ, Δxₙ₋₂)
    b₁₂ = dot(Δxₙ₋₁, Δxₙ₋₂)
    b₁₁ = dot(Δxₙ₋₁, Δxₙ₋₁)
    b₂₂ = dot(Δxₙ₋₂, Δxₙ₋₂)
    den = b₁₁*b₂₂ - b₁₂^2
    μ1 = (b₀₂*b₁₂ - b₀₁*b₂₂) / den
    μ2 = (b₀₁*b₁₂ - b₀₂*b₁₁) / den
    x_gdem = xₙ .+ (Δxₙ .- μ2 .* Δxₙ₋₁)/(1. + μ1 + μ2)
    return x_gdem
end
=#