struct IAPWS95 <: EmpiricHelmholtzModel end

const IAPWS95_params = (
    n00= [-8.3204464837497, 6.6832105275932, 3.00632,0.012436, 0.97315, 1.2795, 0.96956, 0.24873]
    ,gamma00 = [0, 0, 0, 1.28728967, 3.53734222, 7.74073708, 9.24437796,27.5075105]
    
    ,nr1 = [0.12533547935523e-1, 0.78957634722828e1, -0.87803203303561e1,
    0.31802509345418, -0.26145533859358, -0.78199751687981e-2,
    0.88089493102134e-2]
    ,d1 = [1, 1, 1, 2, 2, 3, 4]
    ,t1 = [-0.5, 0.875, 1, 0.5, 0.75, 0.375, 1]
    
    ,nr2 = [-0.66856572307965, 0.20433810950965, -0.66212605039687e-4,
    -0.19232721156002, -0.25709043003438, 0.16074868486251,
    -0.4009282892587e-1, 0.39343422603254e-6, -0.75941377088144e-5,
    0.56250979351888e-3, -0.15608652257135e-4, 0.11537996422951e-8,
    .36582165144204e-6, -.13251180074668e-11, -.62639586912454e-9,
    -0.10793600908932, 0.17611491008752e-1, 0.22132295167546,
    -0.40247669763528, 0.58083399985759, 0.49969146990806e-2,
    -0.31358700712549e-1, -0.74315929710341, 0.47807329915480,
    0.20527940895948e-1, -0.13636435110343, 0.14180634400617e-1,
    0.83326504880713e-2, -0.29052336009585e-1, 0.38615085574206e-1,
    -0.20393486513704e-1, -0.16554050063734e-2, .19955571979541e-2,
    0.15870308324157e-3, -0.16388568342530e-4, 0.43613615723811e-1,
    0.34994005463765e-1, -0.76788197844621e-1, 0.22446277332006e-1,
    -0.62689710414685e-4, -0.55711118565645e-9, -0.19905718354408,
    0.31777497330738, -0.11841182425981]
    ,c2 = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 6, 6,6, 6]
    ,d2 = [1, 1, 1, 2, 2, 3, 4, 4, 5, 7, 9, 10, 11, 13, 15, 1, 2, 2, 2, 3,
    4, 4, 4, 5, 6, 6, 7, 9, 9, 9, 9, 9, 10, 10, 12, 3, 4, 4, 5, 14, 3, 6, 6, 6]
    ,t2= [4, 6, 12, 1, 5, 4, 2, 13, 9, 3, 4, 11, 4, 13, 1, 7, 1, 9, 10, 10, 3, 7, 
    10, 10, 6, 10, 10, 1, 2, 3, 4, 8, 6, 9, 8, 16, 22, 23,23, 10, 50, 44, 46, 50]
    
    ,nr3 = (-0.31306260323435e2, 0.31546140237781e2, -0.25213154341695e4)
    ,t3 = (0, 1, 4)
    ,beta3 = (150, 150, 250)
    ,gamma3 = (1.21, 1.21, 1.25)
    
    ,nr4 =(-0.14874640856724, 0.31806110878444)
    ,b4 =(0.85, 0.95)
    ,C =(28, 32)
    ,D =(700, 800)
)


@inline function _f0(_model::IAPWS95,rho,T)
    model = IAPWS95_params
    delta = rho/322
    tau = 647.096/T
    
    res = log(delta)+model.n00[1]
    res = muladd(model.n00[2],tau,res)
    res = muladd(model.n00[3],log(tau),res)
    
 for i = 4:8
        res = muladd(model.n00[i],log(-expm1(-model.gamma00[i]*tau)),res)
    end
    return res
end

function _fr(_model::IAPWS95,rho,T)
    model = IAPWS95_params
    delta = rho/322.0
    tau = 647.096/T
   
    res=zero(promote_type(typeof(rho),typeof(T)))
    for i = 1:7
        res = muladd((model.nr1[i]*delta^(model.d1[i])),(tau^model.t1[i]),res)
    end

    for i = 1:44
        res = muladd((model.nr2[i]*delta^(model.d2[i])),(tau^model.t2[i]) * exp(-delta^model.c2[i]),res)
    end

    for i = 1:3
           res=res+(model.nr3[i]*delta^3) * (tau^model.t3[i]) * 
           exp(-20*abs2(delta-1)-model.beta3[i]*abs(tau-model.gamma3[i]))
    end
    delta1m2 = (delta-1)^2# 
    tau1m2 = (tau-1)^2
    
    for i = 1:2
        theta = (1-tau) + 0.32*delta1m2^(1/(2*0.3))
        del = theta^2 + 0.2*delta1m2^3.5
        psi = exp(- model.C[i]*delta1m2 - model.D[i]*tau1m2)
        res = res+model.nr4[i]*del^model.b4[i]*psi
    end
    return res
end

@inline _f(model::IAPWS95,rho,T) = _fr(model,rho,T)+_f0(model,rho,T)

const IAPWS_R_corr = 8.3143713575874/R̄

function a_ideal(model::IAPWS95,V,T,z=@SVector [1.0]) 
    Σz = only(z) #single component
    v = V/Σz
    #R value calculated from molecular weight and specific gas constant
     #return 8.3143713575874*T*_f(model, molar_to_weight(1/v,[model.molecularWeight],[1.0]),T)
     #println(molar_to_weight(1/v,[model.molecularWeight],[1.0]))'    
     mass_v =  v*1000.0*0.055508472036052976
     rho = one(mass_v)/mass_v
     return IAPWS_R_corr*_f0(model,rho,T)
end

function a_res(model::IAPWS95,V,T,z=(one(V),)) 
    Σz = only(z) #single component
    v = V/Σz
    #R value calculated from molecular weight and specific gas constant
     #return 8.3143713575874*T*_f(model, molar_to_weight(1/v,[model.molecularWeight],[1.0]),T)
     #println(molar_to_weight(1/v,[model.molecularWeight],[1.0]))'    
     mass_v =  v*1000.0*0.055508472036052976
     rho = one(mass_v)/mass_v
     return IAPWS_R_corr*_fr(model,rho,T)
end

struct IAPWS95Ideal <:IdealModel end

function IAPWS95Ideal(components; verbose=false)
    if only(components) == "water"
        return IAPWS95Ideal()
    else
        return error("IAPWS95 is for water only")
    end
end

function a_ideal(model::IAPWS95Ideal,V,T,z=one(V)) 
    return a_ideal(IAPWS95(),V, T, z)
end

idealmodel(model::IAPWS95) = IAPWS95Ideal()

export IAPWS95,IAPWS95Ideal
