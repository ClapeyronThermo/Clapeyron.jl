struct LJSAFTParam <: EoSParam
    segment::SingleParam{Float64}
    b::PairParam{Float64}
    T_tilde::PairParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
end

abstract type LJSAFTModel <: SAFTModel end
@newmodel LJSAFT LJSAFTModel LJSAFTParam

export LJSAFT
function LJSAFT(components; idealmodel=BasicIdeal, userlocations=String[], ideal_userlocations=String[], verbose=false)
    params,sites = getparams(components, ["SAFT/LJSAFT"]; userlocations=userlocations, verbose=verbose)
    segment = params["m"]

    k = params["k"]
    zeta = params["zeta"]
    
    T_tilde = epsilon_LorentzBerthelot(params["T_tilde"], k)
    params["b"].values .*= 1E-3
    params["b"].values .^= 1/3
    b = sigma_LorentzBerthelot_mod(params["b"],zeta)
    b.values .^= 3
    epsilon_assoc = params["epsilon_assoc"]
    bondvol = params["bondvol"]

    packagedparams = LJSAFTParam(segment, b, T_tilde, epsilon_assoc, bondvol)
    references = ["10.1021/ie9602320"]

    model = LJSAFT(packagedparams, sites, idealmodel; ideal_userlocations=ideal_userlocations, references=references, verbose=verbose)
    return model
end

function a_res(model::LJSAFTModel, V, T, z)
    return @f(a_seg) + @f(a_chain) + @f(a_assoc)
end

function a_seg(model::LJSAFTModel, V, T, z)
    D = LJSAFTconsts.D
    C = LJSAFTconsts.C
    C0 = LJSAFTconsts.C0
    C1 = LJSAFTconsts.C1
    C2 = LJSAFTconsts.C2
    C4 = LJSAFTconsts.C4
  
    Σz = ∑(z)
    m = model.params.segment.values

    T̃ = @f(Tm)
    b̄ = @f(bm)

    Tst = T/T̃
    m̄ = dot(m,z)/Σz
    ρ = Σz/V
    ρst = m̄*b̄*ρ
    η = ρst*π/6*(∑(D[i+3]*Tst^(i/2) for i ∈ -2:1)+D[end]*log(Tst))^3

    A_HS = Tst*(5/3*log(1-η)+(η*(34-33η+4η^2))/(6*(1-η)^2))
    ΔB2 = ∑(C[j+8]*Tst^(j/2) for j ∈ -7:0)
    A0 = ∑(C0[j-1]*ρst^j for j ∈ 2:5)
    A1 = ∑(C1[j-1]*Tst^(-1/2)*ρst^j for j ∈ 2:6)
    A2 = ∑(C2[j-1]*Tst^(-1)*ρst^j for j ∈ 2:6)
    A4 = ∑(C4[j-1]*Tst^(-2)*ρst^j for j ∈ 2:6)

    γ = 1.92907278
    return m̄*(A_HS+exp(-γ*ρst^2)*ρst*Tst*ΔB2+A0+A1+A2+A4)/Tst
end

function Tm(model::LJSAFTModel, V, T, z)
    #x = z/∑(z)
    T̃ = model.params.T_tilde.values
    b = model.params.b.values
    m = model.params.segment.values
    comps = @comps
    return ∑(m[i]*m[j]*z[i]*z[j]*b[i,j]*T̃[i,j] for i ∈ 1:length(z) for j ∈ comps)/∑(m[i]*m[j]*z[i]*z[j]*b[i,j] for i ∈ comps for j ∈ comps)
end

function bm(model::LJSAFTModel, V, T, z)
    comps = @comps
    b = model.params.b.values
    m = model.params.segment.values
    return ∑(m[i]*m[j]*z[i]*z[j]*b[i,j] for i ∈ comps for j ∈ comps)/∑(m[i]*m[j]*z[i]*z[j] for i ∈ comps for j ∈ comps)
end

function a_chain(model::LJSAFTModel, V, T, z)
    m = model.params.segment.values

    return -sum(z[i]*(m[i]-1)*log(@f(g_LJ,i)) for i ∈ @comps)/∑(z)
end

function g_LJ(model::LJSAFTModel, V, T, z, i)

    m = model.params.segment.values
    b = model.params.b.diagvalues
    T̃ = model.params.T_tilde.diagvalues
    Tst = T/T̃[i]
    ρ = 1/V
    ρ̄ = b[i]*m[i]*z[i]*ρ
    a = LJSAFTconsts.a::Matrix{Float64}
    return (1+sum(a[i,j]*ρ̄^i*Tst^(1-j) for i ∈ 1:5 for j ∈ 1:5))
end

function Δ(model::LJSAFTModel, V, T, z, i, j, a, b)
    ∑z = ∑(z)
    m = model.params.segment.values
    _b = model.params.b.values
    T̃ = model.params.T_tilde.values
    Tst = T/T̃[i,j]
    ρ = ∑z/V
    ρ̄ = z[i]*_b[i,j]*m[i]*ρ/∑z
    ϵ_assoc = model.params.epsilon_assoc.values
    κ = model.params.bondvol.values
    b_ = LJSAFTconsts.b
    I = sum(b_[i+1,j+1]*ρ̄^i*Tst^j for i ∈ 0:4 for j ∈ 0:4)/3.84/1e4
    return 4π*(exp(ϵ_assoc[i,j][a,b]/T)-1)*κ[i,j][a,b]*I*_b[i,j]/N_A
end

const LJSAFTconsts = (
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
    D = [0.011117524,-0.076383859,1.080142248, 0.000693129,-0.063920968],
    C = [-0.58544978,0.43102052,0.87361369,-4.13749995,2.90616279,-7.02181962,0.,0.0245987],
    C0 = [2.01546797,-28.17881636,28.28313847,-10.42402873],
    C1 = [-19.58371655,75.62340289,-120.70586598,93.92740328,-27.37737354],
    C2 = [29.34470520,-112.35356937,170.64908980,-123.06669187,34.42288969],
    C4 = [-13.37031968,65.38059570,-115.09233113,88.91973082,-25.62099890]    
   )
