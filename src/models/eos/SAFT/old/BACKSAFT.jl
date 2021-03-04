const BACKSAFT_consts = (
    D1 = [-8.8043,4.164627,-48.203555,140.4362,-195.23339,113.515]
    ,D2 = [2.9396,-6.0865383,40.137956,-76.230797,-133.70055,860.25349,-1535.3224,1221.4261,-409.10539]
    ,D3 = [-2.8225,4.7600148,11.257177,-66.382743,69.248785]
    ,D4 = [0.34,-3.1875014,12.231796,-12.110681]
)

function a_res(model::BACKSAFTFamily, z, v, T)
    ahcb = a_hcb(model,z,v,T)
    adisp = a_disp(model,z,v,T)
    achain = a_chain(model,z,v,T)
    return  ahcb + achain + (1.75*(achain/ahcb)+1)*adisp
end

function a_hcb(model::BACKSAFTFamily, z, v, T)
    component = model.components[1]
    α = model.params.alpha[component]
    m = model.params.segment[component]
    η = ζn(model,z,v,T, 3)
    return m*(α^2/(1-η)^2-(α^2-3α)/(1-η)-(1-α^2)*log(1-η)-3α)
end

function a_disp(model::BACKSAFTFamily, z, v, T)
    component = model.components[1]
    m = model.params.segment[component]
    c = model.params.c[component]
    u = model.params.epsilon[component]*(1+c/T)
    η = ζn(model,z,v,T,3)
    τ = 0.74048
    D1 = BACKSAFT_consts.D1
    D2 = BACKSAFT_consts.D2
    D3 = BACKSAFT_consts.D3
    D4 = BACKSAFT_consts.D4
    A1 = sum(D1[j]*(u/T)*(η/τ)^j for j in 1:6)
    A2 = sum(D2[j]*(u/T)^2*(η/τ)^j for j in 1:9)
    A3 = sum(D3[j]*(u/T)^3*(η/τ)^j for j in 1:5)
    A4 = sum(D4[j]*(u/T)^4*(η/τ)^j for j in 1:4)
    return m*(A1+A2+A3+A4)
end

function d(model::BACKSAFTFamily, z, v, T, component)
    ϵ = model.params.epsilon[component]
    σ = model.params.sigma[component]
    return σ * (1 - 0.12exp(-3ϵ/T))
end

function ζn(model::BACKSAFTFamily, z, v, T, n)
    x = z/sum(z[i] for i in model.components)
    m = model.params.segment
    return N_A*sum(z[i] for i in model.components)*π/6/v * sum(x[i]*m[i]*d(model,z,v,T, i)^n for i in model.components)
end

function a_chain(model::BACKSAFTFamily, z, v, T)
    component = model.components[1]
    m = model.params.segment[component]
    return (1-m)*log(g_hcb(model,z,v,T))
end

function g_hcb(model::BACKSAFTFamily, z, v, T)
    component = model.components[1]
    α = model.params.alpha[component]
    η = ζn(model,z,v,T, 3)
    return 1/(1-η)+3*(1+α)*α*η/((1-η)^2*(1+3α))+3*η^2*α^2/((1-η)^3*(1+3α))
end
