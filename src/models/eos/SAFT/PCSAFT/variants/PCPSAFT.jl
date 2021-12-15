struct PCPSAFTParam <: EoSParam
    Mw::SingleParam{Float64}
    segment::SingleParam{Float64}
    mu::SingleParam{Float64}
    sigma::PairParam{Float64}
    epsilon::PairParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
end

abstract type PCPSAFTModel <: PCSAFTModel end
@newmodel PCPSAFT PCPSAFTModel PCPSAFTParam

export PCPSAFT
function PCPSAFT(components; idealmodel=BasicIdeal, userlocations=String[], ideal_userlocations=String[], verbose=false)
    params,sites = getparams(components, [ "SAFT/PCSAFT/PCPSAFT"]; userlocations=userlocations, verbose=verbose)
    
    segment = params["m"]
    k = params["k"]
    Mw = params["Mw"]
    params["mu"].values .*=3.33564e-30/√(4*π*ϵ_0)
    mu = params["mu"]
    params["sigma"].values .*= 1E-10
    sigma = sigma_LorentzBerthelot(params["sigma"])
    epsilon = epsilon_LorentzBerthelot(params["epsilon"], k)
    epsilon_assoc = params["epsilon_assoc"]
    bondvol = params["bondvol"]
    packagedparams = PCPSAFTParam(Mw, segment, mu, sigma, epsilon, epsilon_assoc, bondvol)
    references = ["10.1021/ie020753p"]

    model = PCPSAFT(packagedparams, sites, idealmodel; ideal_userlocations=ideal_userlocations, references=references, verbose=verbose)
    return model
end

function a_res(model::PCPSAFTModel, V, T, z)
    return @f(a_hc) + @f(a_disp) + @f(a_assoc) + @f(a_multipole)
end

function a_multipole(model::PCPSAFTModel, V, T, z)
    if sum(model.params.mu.values)==0.
        return 0.
    else
        return @f(a_2)/(1-@f(a_3)/@f(a_2))
    end
end

function a_2(model::PCPSAFTModel, V, T, z)
    μ = model.params.mu.values
    ϵ = model.params.epsilon.diagvalues
    σ = model.params.sigma.diagvalues
    m = model.params.segment.values

    μ̄ = @. μ/√(m*σ^3*ϵ*k_B)
    T̄ = @. T/ϵ

    σ = model.params.sigma.values

    ρ = N_A*∑(z)/V
    x = z/∑(z)
    return -π*ρ*∑(x[i]*x[j]*μ̄[i]^2*μ̄[j]^2*σ[i,i]^3*σ[j,j]^3/(σ[i,j]^3)/T̄[i]/T̄[j]*@f(J₂,i,j) for i ∈ @comps, j ∈ @comps)
end

function a_3(model::PCPSAFTModel, V, T, z)
    μ = model.params.mu.values
    ϵ = model.params.epsilon.diagvalues
    σ = model.params.sigma.diagvalues
    m = model.params.segment.values

    μ̄ = @. μ/√(m*σ^3*ϵ*k_B)
    T̄ = @. T/ϵ

    σ = model.params.sigma.values

    ρ = N_A*∑(z)/V
    x = z/∑(z)
    return -4/3*π^2*ρ^2*∑(x[i]*x[j]*x[k]*μ̄[i]^2*μ̄[j]^2*μ̄[k]^2*σ[i,i]^3*σ[j,j]^3*σ[k,k]^3/(σ[i,j]*σ[i,k]*σ[j,k])/T̄[i]/T̄[j]/T̄[k]*@f(J₃,i,j,k) for i ∈ @comps, j ∈ @comps, k ∈ @comps)
end

function J₂(model::PCPSAFTModel, V, T, z, i, j)
    η = @f(ζ,3)

    a = PCPSAFTconsts.corr1
    b = PCPSAFTconsts.corr2
    
    m = model.params.segment.values
    m̄ = min(2,√(m[i]*m[j]))

    ϵ = model.params.epsilon.values[i,j]
    T̄ = T/ϵ

    return ∑(((a[n+1,1] + (m̄-1)/m̄*a[n+1,2] + (m̄-1)/m̄*(m̄-2)/m̄*a[n+1,3])+(b[n+1,1] + (m̄-1)/m̄*b[n+1,2] + (m̄-1)/m̄*(m̄-2)/m̄*b[n+1,3])/T̄) * η^n for n = 0:4)
end

function J₃(model::PCPSAFTModel, V, T, z, i, j, k)
    η = @f(ζ,3)

    c = PCPSAFTconsts.corr3
    
    m = model.params.segment.values
    m̄ = min(2,(m[i]*m[j]*m[k])^(1/3))

    return ∑((c[n+1,1] + (m̄-1)/m̄*c[n+1,2] + (m̄-1)/m̄*(m̄-2)/m̄*c[n+1,3])*η^n for n = 0:4)
end

const PCPSAFTconsts = (
    corr1 =
    [0.3043504	0.9534641	-1.161008;
    -0.1358588	-1.8396383	4.5258607;
    1.4493329	2.013118	0.9751222;
    0.3556977	-7.3724958	-12.281038;
    -2.0653308	8.2374135	5.9397575],

    corr2 =
    [0.2187939	-0.5873164	3.4869576;
    -1.1896431	1.2489132	-14.915974;
    1.1626889	-0.508528	15.372022;
    0	0	0;
    0	0	0],

    corr3 =
    [-0.0646774	-0.9520876	-0.6260979;
    0.1975882	2.9924258	1.2924686;
    -0.8087562	-2.3802636	1.6542783;
    0.6902849	-0.2701261	-3.4396744;
    0	0	0]
)