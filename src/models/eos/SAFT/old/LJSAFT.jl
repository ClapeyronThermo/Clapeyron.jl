function a_res(model::LJSAFTFamily, z,Vol,Temp)
    return a_LJ(model,z,Vol,Temp)+a_chain(model,z,Vol,Temp)+a_assoc(model,z,Vol,Temp)
end

function a_LJ(model::LJSAFTFamily, z,Vol,Temp)
    D   = [0.011117524,-0.076383859,1.080142248, 0.000693129,-0.063920968]
    C   = [-0.58544978,0.43102052,0.87361369,-4.13749995,2.90616279,-7.02181962,0.,0.0245987]
    C0  = [2.01546797,-28.17881636,28.28313847,-10.42402873]
    C1  = [-19.58371655,75.62340289,-120.70586598,93.92740328,-27.37737354]
    C2  = [29.34470520,-112.35356937,170.64908980,-123.06669187,34.42288969]
    C4  = [-13.37031968,65.38059570,-115.09233113,88.91973082,-25.62099890]
    x   = z/sum(z[i] for i in model.components)
    m   = model.params.segment

    T̃    = Tm(model,z,Vol,Temp)
    b̄    = bm(model,z,Vol,Temp)

    T̄    = Temp/T̃
    m̄    = sum(x[i]*m[i] for i in model.components)
    ρ    = sum(z)/Vol
    ρ̄    = m̄*b̄*ρ
    η    = ρ̄*π/6*(sum(D[i+3]*T̄^(i/2) for i in -2:1)+D[end]*log(T̄))^3

    A_HS = T̄*(5/3*log(1-η)+(η*(34-33η+4η^2))/(6*(1-η)^2))
    ΔB2  = sum(C[j+8]*T̄^(j/2) for j in -7:0)
    A0   = sum(C0[j-1]*ρ̄^j for j in 2:5)
    A1   = sum(C1[j-1]*T̄^(-1/2)*ρ̄^j for j in 2:6)
    A2   = sum(C2[j-1]*T̄^(-1)*ρ̄^j for j in 2:6)
    A4   = sum(C4[j-1]*T̄^(-2)*ρ̄^j for j in 2:6)

    γ    = 1.92907278
    return m̄*(A_HS+exp(-γ*ρ̄^2)*ρ̄*T̄*ΔB2+A0+A1+A2+A4)/T̄
end

function Tm(model::LJSAFTFamily,z,Vol,Temp)
    components = model.components
    x   = z/sum(z[i] for i in model.components)
    T̃    = model.params.T
    b    = model.params.b
    m    = model.params.segment
    return sum(m[i]*m[j]*x[i]*x[j]*b[union(i,j)]*T̃[union(i,j)] for i in components for j in components)/sum(m[i]*m[j]*x[i]*x[j]*b[union(i,j)] for i in components for j in components)
end

function bm(model::LJSAFTFamily,z,Vol,Temp)
    components = model.components
    x   = z/sum(z[i] for i in model.components)
    b    = model.params.b
    m    = model.params.segment
    return sum(m[i]*m[j]*x[i]*x[j]*b[union(i,j)] for i in components for j in components)/sum(m[i]*m[j]*x[i]*x[j] for i in components for j in components)
end

function a_chain(model::LJSAFTFamily, z,Vol,Temp)
    x       = z/sum(z[i] for i in model.components)
    m̄      = sum(model.params.segment[i]*x[i] for i in model.components)
    return -log(g_LJ(model,z,Vol,Temp))*(m̄-1)
end

function g_LJ(model::LJSAFTFamily,z,Vol,Temp)
    x   = z/sum(z[i] for i in model.components)
    m   = model.params.segment

    T̃    = Tm(model,z,Vol,Temp)
    b̄    = bm(model,z,Vol,Temp)

    T̄    = Temp/T̃
    m̄    = sum(x[i]*m[i] for i in model.components)
    ρ    = sum(z)/Vol
    ρ̄    = m̄*b̄*ρ
    a    = [0.49304346593882 2.1528349894745 -15.955682329017 24.035999666294 -8.6437958513990;
           -0.47031983115362 1.1471647487376  37.889828024211 -84.667121491179 39.643914108411;
            5.0325486243620 -25.915399226419 -18.862251310090 107.63707381726 -66.602649735720;
           -7.3633150434385  51.553565337453 -40.519369256098 -38.796692647218 44.605139198378;
            2.9043607296043 -24.478812869291  31.500186765040 -5.3368920371407 -9.5183440180133]
    return (1+sum(a[i,j]*ρ̄^i*T̄^(1-j) for i in 1:5 for j in 1:5))
end

function a_assoc(model::LJSAFTFamily, z, v, T)
    x = z/sum(z[i] for i in model.components)
    n_sites = model.params.n_sites
    X_iA = X_assoc(model,z,v,T)
    return sum(x[i]*sum(n_sites[i][a]*(log(X_iA[i,a])+(1-X_iA[i,a])/2) for a in keys(model.params.n_sites[i])) for i in model.components)
end

function X_assoc(model::LJSAFTFamily, z, v, T)
    x = z/sum(z[i] for i in model.components)
    ρ = sum(z[i] for i in model.components)/v
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

function Δ(model::LJSAFTFamily, z, Vol, Temp, i, j, a, b)
    x   = z/sum(z[i] for i in model.components)
    m   = model.params.segment

    T̃    = Tm(model,z,Vol,Temp)
    b̄    = bm(model,z,Vol,Temp)

    T̄    = Temp/T̃
    m̄    = sum(x[i]*m[i] for i in model.components)
    ρ    = sum(z)/Vol
    ρ̄    = m̄*b̄*ρ
    ϵ_assoc = model.params.epsilon_assoc[Set([(i,a),(j,b)])]
    κ = model.params.bond_vol[Set([(i,a),(j,b)])]

    b    = [-0.03915181 0.08450471 0.06889053 -0.01034279  0.5728662e-3;
             -0.5915018  0.9838141 -0.4862279   0.1029708 -0.6919154e-2;
               1.908368  -3.415721   2.124052  -0.4298159  0.02798384;
             -0.7957312  0.7187330 -0.9678804   0.2431675 -0.01644710;
             -0.9399577   2.314054 -0.4877045  0.03932058 -0.1600850e-2]

    I = sum(b[i+1,j+1]*ρ̄^i*T̄^j for i in 0:4 for j in 0:4)/3.84/1e4
    return 4π*(exp(ϵ_assoc/Temp)-1)*κ*I*b̄
end
