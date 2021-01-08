function a_res(model::softSAFTFamily, z,Vol,Temp)
    return a_LJ(model,z,Vol,Temp)+a_chain(model,z,Vol,Temp)+a_assoc(model,z,Vol,Temp)
end

function a_LJ(model::softSAFTFamily, z,Vol,Temp)
    ϵ̄    = ϵm(model,z,Vol,Temp)
    σ̄    = σm(model,z,Vol,Temp)
    m̄    = sum(z[i]*model.params.segment[i] for i in model.components)/sum(z)
    T̄    = Temp/ϵ̄
    ρ̄    = ρs(model,z,Vol,Temp)*σ̄^3
    γ    = 3
    F    = exp(-γ*ρ̄^2)
    x    = [ 0.8623085097507421, 2.976218765822098,-8.402230115796038,0.1054136629203555,
            -0.8564583828174598, 1.582759470107601,0.7639421948305453, 1.753173414312048,
            2.798291772190376e3,-4.8394220260857657e-2,0.9963265197721935,-3.698000291272493e1,
            2.084012299434647e1, 8.305402124717285e1,-9.574799715203068e2,-1.477746229234994e2,
            6.398607852471505e1, 1.603993673294834e1, 6.805916615864377e1,-2.791293578795945e3,
            -6.245128304568454, -8.116836104958410e3, 1.488735559561229e1,-1.059346754655084e4,
            -1.131607632802822e2,-8.867771540418822e3,-3.986982844450543e1,-4.689270299917261e3,
            2.593535277438717e2,-2.694523589434903e3,-7.218487631550215e2, 1.721802063863269e2]
    a    = []
    b    = []
    G    = []

    append!(a,x[1]*T̄+x[2]*√(T̄)+x[3]+x[4]/T̄+x[5]/T̄^2)
    append!(a,x[6]*T̄+x[7]+x[8]/T̄+x[9]/T̄^2)
    append!(a,x[11]+x[10]*T̄+x[12]/T̄)
    append!(a,x[13])
    append!(a,x[14]/T̄+x[15]/T̄^2)
    append!(a,x[16]/T̄)
    append!(a,x[17]/T̄+x[18]/T̄^2)
    append!(a,x[19]/T̄^2)

    append!(b,x[20]/T̄^2+x[21]/T̄^3)
    append!(b,x[22]/T̄^2+x[23]/T̄^4)
    append!(b,x[24]/T̄^2+x[25]/T̄^3)
    append!(b,x[26]/T̄^2+x[27]/T̄^4)
    append!(b,x[28]/T̄^2+x[29]/T̄^3)
    append!(b,x[30]/T̄^2+x[31]/T̄^3+x[32]/T̄^4)

    append!(G,(1-F)/(2γ))
    append!(G,-(F*ρ̄^2-2G[1])/(2γ))
    append!(G,-(F*ρ̄^4-4G[2])/(2γ))
    append!(G,-(F*ρ̄^6-6G[3])/(2γ))
    append!(G,-(F*ρ̄^8-8G[4])/(2γ))
    append!(G,-(F*ρ̄^10-10G[5])/(2γ))
    return m̄*(sum(a[i]*ρ̄^i/i for i in 1:8)+sum(b[i]*G[i] for i in 1:6))/T̄
end

function ϵm(model::softSAFTFamily,z,Vol,Temp)
    components = model.components
    x   = z/sum(z[i] for i in model.components)
    ϵ    = model.params.epsilon
    σ    = model.params.sigma
    m    = model.params.segment
    return sum(m[i]*m[j]*x[i]*x[j]*σ[union(i,j)]^3*ϵ[union(i,j)] for i in components for j in components)/sum(m[i]*m[j]*x[i]*x[j]*σ[union(i,j)]^3 for i in components for j in components)
end

function σm(model::softSAFTFamily,z,Vol,Temp)
    components = model.components
    x   = z/sum(z[i] for i in model.components)
    σ    = model.params.sigma
    m    = model.params.segment
    return (sum(m[i]*m[j]*x[i]*x[j]*σ[union(i,j)]^3 for i in components for j in components)/sum(m[i]*m[j]*x[i]*x[j] for i in components for j in components))^(1/3)
end

function ρs(model::softSAFTFamily, z,Vol,Temp)
    NParticles = N_A*sum(z[i] for i in model.components)
    x   = z/sum(z[i] for i in model.components)
    m   = model.params.segment
    m̄   = sum(x[i]*m[i] for i in model.components)
    return NParticles/Vol*m̄
end

function xS(model::softSAFTFamily, z,Vol,Temp,i)
    x   = z/sum(z[i] for i in model.components)
    m   = model.params.segment
    m̄   = sum(x[j]*m[j] for j in model.components)
    return x[i]*m[i]/m̄
end

function a_chain(model::softSAFTFamily, z,Vol,Temp)
    x       = z/sum(z[i] for i in model.components)
    m̄      = sum(model.params.segment[i]*x[i] for i in model.components)
    return -log(y_LJ(model,z,Vol,Temp))*(m̄-1)
end

function y_LJ(model::softSAFTFamily,z,Vol,Temp)
    ϵ̄    = ϵm(model,z,Vol,Temp)
    gLJ  = g_LJ(model,z,Vol,Temp)
    return gLJ*exp(-ϵ̄/Temp)
end

function g_LJ(model::softSAFTFamily,z,Vol,Temp)
    ϵ̄    = ϵm(model,z,Vol,Temp)
    σ̄    = σm(model,z,Vol,Temp)
    T̄    = Temp/ϵ̄
    ρ̄    = ρs(model,z,Vol,Temp)*σ̄^3
    a    = [0.49304346593882 2.1528349894745 -15.955682329017 24.035999666294 -8.6437958513990;
           -0.47031983115362 1.1471647487376  37.889828024211 -84.667121491179 39.643914108411;
            5.0325486243620 -25.915399226419 -18.862251310090 107.63707381726 -66.602649735720;
           -7.3633150434385  51.553565337453 -40.519369256098 -38.796692647218 44.605139198378;
            2.9043607296043 -24.478812869291  31.500186765040 -5.3368920371407 -9.5183440180133]
    return 1+sum(a[i,j]*ρ̄^i*T̄^(1-j) for i in 1:5 for j in 1:5)
end

function a_assoc(model::softSAFTFamily, z, v, T)
    x = z/sum(z[i] for i in model.components)
    n_sites = model.params.n_sites
    X_iA = X_assoc(model,z,v,T)
    return sum(x[i]*sum(n_sites[i][a]*(log(X_iA[i,a])+(1-X_iA[i,a])/2) for a in keys(model.params.n_sites[i])) for i in model.components)
end

function X_assoc(model::softSAFTFamily, z, v, T)
    x = z/sum(z[i] for i in model.components)
    ρ = N_A*sum(z[i] for i in model.components)/v
    n_sites = model.params.n_sites
    X_iA = Dict()
    X_iA_old = Dict()
    tol = 1.
    iter = 1

    while tol > 1e-12 && iter < 100
        for i in model.components
            for a in keys(model.params.n_sites[i])
                A = 0.
                for j in model.components
                    B = 0
                    for b in keys(model.params.n_sites[j])
                        if haskey(model.params.epsilon_assoc,Set([(i,a),(j,b)]))
                            if iter!=1
                                B+=n_sites[j][b]*X_iA_old[j,b]*Δ(model,z,v,T,i,j,a,b)
                            else
                                B+=n_sites[j][b]*Δ(model,z,v,T,i,j,a,b)
                            end
                        end
                    end
                    A += ρ*x[j]*B
                end
                if iter == 1
                    X_iA[i,a] =0.5+0.5*(1+A)^-1
                else
                    X_iA[i,a] =0.5*X_iA_old[i,a]+0.5*(1+A)^-1
                end
            end
        end
        if iter == 1
            tol = sqrt(sum(sum((1. -X_iA[i,a])^2 for a in keys(model.params.n_sites[i])) for i in model.components))
        else
            tol = sqrt(sum(sum((X_iA_old[i,a] -X_iA[i,a])^2 for a in keys(model.params.n_sites[i])) for i in model.components))
        end
        X_iA_old = deepcopy(X_iA)
        iter += 1
    end

    return X_iA
end

function Δ(model::softSAFTFamily, z, Vol, Temp, i, j, a, b)
    ϵ̄    = ϵm(model,z,Vol,Temp)
    σ̄    = σm(model,z,Vol,Temp)
    T̄    = Temp/ϵ̄
    ρ̄    = ρs(model,z,Vol,Temp)*σ̄^3
    ϵ_assoc = model.params.epsilon_assoc[Set([(i,a),(j,b)])]
    κ = model.params.bond_vol[Set([(i,a),(j,b)])]

    b    = [-0.03915181 0.08450471 0.06889053 -0.01034279  0.5728662e-3;
             -0.5915018  0.9838141 -0.4862279   0.1029708 -0.6919154e-2;
               1.908368  -3.415721   2.124052  -0.4298159  0.02798384;
             -0.7957312  0.7187330 -0.9678804   0.2431675 -0.01644710;
             -0.9399577   2.314054 -0.4877045  0.03932058 -0.1600850e-2]

    I = sum(b[i+1,j+1]*ρ̄^i*T̄^j for i in 0:4 for j in 0:4)/3.84/1e4
    return 4π*(exp(ϵ_assoc/Temp)-1)*κ*I
end
