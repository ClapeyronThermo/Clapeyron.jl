struct softSAFTParam <: EoSParam
    segment::SingleParam{Float64}
    sigma::PairParam{Float64}
    epsilon::PairParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
end

abstract type softSAFTModel <: SAFTModel end
@newmodel softSAFT softSAFTModel softSAFTParam

export softSAFT
function softSAFT(components; idealmodel=BasicIdeal, userlocations=String[], ideal_userlocations=String[], verbose=false)
    params,sites = getparams(components, ["SAFT/softSAFT"]; userlocations=userlocations, verbose=verbose)
    
    segment = params["m"]
    k = params["k"]
    params["sigma"].values .*= 1E-10
    sigma = sigma_LorentzBerthelot(params["sigma"])
    epsilon = epsilon_LorentzBerthelot(params["epsilon"], k)
    epsilon_assoc = params["epsilon_assoc"]
    bondvol = params["bondvol"]

    packagedparams = softSAFTParam(segment, sigma, epsilon, epsilon_assoc, bondvol)
    references = ["todo"]

    model = softSAFT(packagedparams, sites, idealmodel; ideal_userlocations=ideal_userlocations, references=references, verbose=verbose)
    return model
end

function a_res(model::softSAFTModel, V, T, z)
    return @f(a_LJ) + @f(a_chain) + @f(a_assoc)
end

function a_LJ(model::softSAFTModel, V, T, z)
    m = model.params.segment.values
    Σz = sum(z)
    m̄ = dot(z,m)/Σz
    ϵ̄ = @f(ϵ_m)
    σ̄ = @f(σ_m)
    ρ̄ = @f(ρ_S)*σ̄^3
    T̄ = T/ϵ̄
    γ = 3
    F = exp(-γ*ρ̄^2)
    x = softSAFTconsts.x
    T2 = T̄*T̄
    T3 = T2*T̄
    T4 = T2*T2
    a = (
         x[1]*T̄+x[2]*√(T̄)+x[3]+x[4]/T̄+x[5]/T̄^2,
         x[6]*T̄+x[7]+x[8]/T̄+x[9]/T̄^2,
         x[11]+x[10]*T̄+x[12]/T̄,
         x[13],
         x[14]/T̄+x[15]/T̄^2,
         x[16]/T̄,
         x[17]/T̄+x[18]/T̄^2,
         x[19]/T̄^2,
        )
    b = (
         x[20]/T2+x[21]/T3,
         x[22]/T2+x[23]/T4,
         x[24]/T2+x[25]/T3,
         x[26]/T2+x[27]/T4,
         x[28]/T2+x[29]/T3,
         x[30]/T2+x[31]/T3+x[32]/T4,
        )
    G1 = (1-F)/(2γ)
    G2 = -(F*ρ̄ ^2 - 2*G1) / 2γ
    G3 = -(F*ρ̄ ^4 - 4*G2) / 2γ
    G4 = -(F*ρ̄ ^6 - 6*G3) / 2γ
    G5 = -(F*ρ̄ ^8 - 8*G4) / 2γ
    G6 = -(F*ρ̄ ^10 - 10*G5) / 2γ
    G = (G1,G2,G3,G4,G5,G6)
  
    return m̄*(∑(a[i]*ρ̄ ^i/i for i ∈ 1:8)+∑(b[i]*G[i] for i ∈ 1:6))/T̄
end
function ϵ_m(model::softSAFTModel, V, T, z)
    comps = @comps
    ϵ = model.params.epsilon.values
    σ = model.params.sigma.values
    m = model.params.segment.values
    return sum(m[i]*m[j]*z[i]*z[j]*σ[i,j]^3*ϵ[i,j] for i ∈ comps for j ∈ comps)/sum(m[i]*m[j]*z[i]*z[j]*σ[i,j]^3 for i ∈ comps for j ∈ comps)
end

function σ_m(model::softSAFTModel, V, T, z)
    comps = @comps
    σ = model.params.sigma.values
    m = model.params.segment.values
    return (sum(m[i]*m[j]*z[i]*z[j]*σ[i,j]^3 for i ∈ comps for j ∈ comps)/sum(m[i]*m[j]*z[i]*z[j] for i ∈ comps for j ∈ comps))^(1/3)
end

function ρ_S(model::softSAFTModel, V, T, z)
    ∑z = ∑(z)
    N = N_A*∑z
    m = model.params.segment.values
    m̄ = dot(z,m)/∑z
    return N/V*m̄
end

function a_chain(model::softSAFTModel, V, T, z)
    m = model.params.segment.values
    m̄ = dot(z,m)/∑(z)
    return -log(@f(y_LJ))*(m̄-1)
end

function y_LJ(model::softSAFTModel, V, T, z)
    ϵ̄ = @f(ϵ_m)
    gLJ = @f(g_LJ)
    return gLJ*exp(-ϵ̄/T)
end

function g_LJ(model::softSAFTModel, V, T, z)
    ϵ̄ = @f(ϵ_m)
    σ̄ = @f(σ_m)
    T̄ = T/ϵ̄
    ρ̄ = @f(ρ_S)*σ̄^3
    a = LJSAFTconsts.a
    return 1+sum(a[i,j]*ρ̄^i*T̄^(1-j) for i ∈ 1:5 for j ∈ 1:5)
end

#=
function X(model::softSAFTModel, V, T, z)
    _1 = one(V+T+first(z))
    ∑z = ∑(z)
    x = z/∑z
    ρ = N_A*∑z/V
    n = model.sites.n_sites
    itermax = 500
    dampingfactor = 0.5
    error = 1.
    tol = model.absolutetolerance
    iter = 1
    X_ = [[_1 for a ∈ @sites(i)] for i ∈ @comps]
    X_old = deepcopy(X_)
    while error > tol
        iter > itermax && error("X has failed to converge after $itermax iterations")
        for i ∈ @comps, a ∈ @sites(i)
            rhs = 1/(1+∑(ρ*x[j]*∑(n[j][b]*X_old[j][b]*@f(Δ,i,j,a,b) for b ∈ @sites(j)) for j ∈ @comps))
            X_[i][a] = (1-dampingfactor)*X_old[i][a] + dampingfactor*rhs
        end
        error = sqrt(∑(∑((X_[i][a] - X_old[i][a])^2 for a ∈ @sites(i)) for i ∈ @comps))
        for i = 1:length(X_)
            X_old[i] .= X_[i]
        end
        iter += 1
    end
    return X_
end
=#
function Δ(model::softSAFTModel, V, T, z, i, j, a, b)
    ϵ̄ = @f(ϵ_m)
    σ̄ = @f(σ_m)
    T̄ = T/ϵ̄
    ρ̄ = @f(ρ_S)*σ̄^3
    ϵ_assoc = model.params.epsilon_assoc.values
    κ = model.params.bondvol.values
    b_ = LJSAFTconsts.b

    I = sum(b_[i+1,j+1]*ρ̄^i*T̄^j for i ∈ 0:4 for j ∈ 0:4)/3.84/1e4
    return 4π*(exp(ϵ_assoc[i,j][a,b]/T)-1)*κ[i,j][a,b]*I
end

const softSAFTconsts =
(
 x = [ 0.8623085097507421, 2.976218765822098,-8.402230115796038,0.1054136629203555,
      -0.8564583828174598, 1.582759470107601,0.7639421948305453, 1.753173414312048,
      2.798291772190376e3,-4.8394220260857657e-2,0.9963265197721935,-3.698000291272493e1,
      2.084012299434647e1, 8.305402124717285e1,-9.574799715203068e2,-1.477746229234994e2,
      6.398607852471505e1, 1.603993673294834e1, 6.805916615864377e1,-2.791293578795945e3,
      -6.245128304568454, -8.116836104958410e3, 1.488735559561229e1,-1.059346754655084e4,
      -1.131607632802822e2,-8.867771540418822e3,-3.986982844450543e1,-4.689270299917261e3,
      2.593535277438717e2,-2.694523589434903e3,-7.218487631550215e2, 1.721802063863269e2],

 a = [0.49304346593882 2.1528349894745 -15.955682329017 24.035999666294 -8.6437958513990;
      -0.47031983115362 1.1471647487376  37.889828024211 -84.667121491179 39.643914108411;
      5.0325486243620 -25.915399226419 -18.862251310090 107.63707381726 -66.602649735720;
      -7.3633150434385  51.553565337453 -40.519369256098 -38.796692647218 44.605139198378;
      2.9043607296043 -24.478812869291  31.500186765040 -5.3368920371407 -9.5183440180133],

 b = [-0.03915181 0.08450471 0.06889053 -0.01034279  0.5728662e-3;
      -0.5915018  0.9838141 -0.4862279   0.1029708 -0.6919154e-2;
      1.908368  -3.415721   2.124052  -0.4298159  0.02798384;
      -0.7957312  0.7187330 -0.9678804   0.2431675 -0.01644710;
      -0.9399577   2.314054 -0.4877045  0.03932058 -0.1600850e-2],
)