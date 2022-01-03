function x0_volume_liquid(model::SAFTgammaMieModel,T,z)
    v_lb = lb_volume(model,z)
    return v_lb*2.0
end
function lb_volume(model::SAFTgammaMieModel, z = SA[1.0])
    vk  = model.groups.n_flattenedgroups
    seg = model.params.segment.values
    S   = model.params.shapefactor.values
    σᵢᵢ = model.params.sigma.diagvalues
    val = π/6*N_A*sum(z[i]*sum(vk[i][k]*seg[k]*S[k]*σᵢᵢ[k]^3 for k in @groups(i)) for i in @comps)
    return val
end

function T_scale(model::SAFTgammaMieModel,z=SA[1.0])
    return T_scale(model.vrmodel,z)
end

function p_scale(model::SAFTgammaMieModel,z=SA[1.0])
    V = zero(first(z))
    T = zero(first(z))
    σ̄3 = @f(σ3x)
    ϵ̄ = T_scale(model,z)
    val    = σ̄3*N_A/R̄/ϵ̄
    return 1/val
end

function a_res(model::SAFTgammaMieModel, V, T, z)
    _data = @f(data)
    return @f(a_hs,_data) + @f(a_disp,_data) + @f(a_chain,_data) + @f(a_assoc,_data)
end

function a_mono(model::SAFTgammaMieModel, V, T, z,_data = @f(data))
    return @f(a_hs,_data) + @f(a_disp,_data)
end

function data(model::SAFTgammaMieModel, V, T, z)
    _d_gc = @f(d_gc)
    _d = @f(d,_d_gc)
    ζi = @f(ζ0123,_d_gc)
    _ζ_X,σ3x = @f(ζ_X_σ3,_d_gc)
    _ρ_S = @f(ρ_S)
    _ζst = _ζst = σ3x*_ρ_S*π/6
    return (_d_gc,(_d,_ρ_S,ζi,_ζ_X,_ζst,σ3x))
end

function d_gc(model::SAFTgammaMieModel, V, T, z)
    _ϵ = model.params.epsilon.diagvalues
    _σ = model.params.sigma.diagvalues
    _λa = model.params.lambda_a.diagvalues
    _λr = model.params.lambda_r.diagvalues
    u = SAFTγMieconsts.u
    w = SAFTγMieconsts.w
    _0 = zero(T+first(z))
    n = length(_ϵ)
    _d_gc = zeros(typeof(_0),n)
    for k ∈ 1:n
        ϵ = _ϵ[k]
        σ = _σ[k]
        λa = _λa[k]
        λr = _λr[k]
        θ = Cλ(model.vrmodel,V,T,z,λa,λr)*ϵ/T
        di = _0
        λrinv = 1/λr
        λaλr = λa/λr
        for j ∈ 1:5
            di += w[j]*(θ/(θ+u[j]))^λrinv*(exp(θ*(1/(θ/(θ+u[j]))^λaλr-1))/(u[j]+θ)/λr)
        end
        _d_gc[k] = σ*(1-di)
    end
    return _d_gc
end

function d(model::SAFTgammaMieModel,V,T,z,_d_gc = @f(d_gc))
    _d = zeros(eltype(_d_gc),length(z))
    _0 = zero(eltype(_d))
    _zz = model.mie_zfractions.z.values
    for i ∈ @comps
        _z = _zz[i]
        m = length(_z)
        di = _0
        for k ∈ 1:m
            zk = _z[k]
            iszero(zk) && continue
            dk = _d_gc[k]
            di += zk*zk*dk^3
            for l ∈ 1:k - 1
                di += 0.25*zk*_z[l]*(dk + _d_gc[l])^3
            end
        end
        _d[i] = cbrt(di)
    end
    #d̂[i] = cbrt(∑(∑(@f(ẑ,i,k)*@f(ẑ,i,l)*@f(d,k,l)^3 for l ∈ @groups) for k ∈ @groups))
    return _d
end

function ζ0123(model::SAFTgammaMieModel, V, T, z,_d_gc=@f(d_gc))
    vrmodel = model.vrmodel
    m = vrmodel.params.segment.values
    m̄ = dot(z, m)
    m̄inv = 1/m̄
    _0 = zero(V+T+first(z))
    ζ0,ζ1,ζ2,ζ3 = _0,_0,_0,_0
    _ϵ = model.params.epsilon.diagvalues
    σ = model.params.sigma.diagvalues
    _λr = model.params.lambda_r.diagvalues
    _λa = model.params.lambda_a.diagvalues
    u = SAFTVRMieconsts.u
    _mi = model.params.mixedsegment.values
    w = SAFTVRMieconsts.w
    for i ∈ @groups
        λa = _λa[i]
        λr = _λr[i]
        ϵ = _ϵ[i]  
        di =_d_gc[i]
        mi = _mi[i]
        xS = dot(mi,z)*m̄inv
        ζ0 += xS
        ζ1 += xS*di
        ζ2 += xS*di*di
        ζ3 += xS*di*di*di
    end
    c = π/6*N_A*m̄/V
    ζ0,ζ1,ζ2,ζ3 = c*ζ0,c*ζ1,c*ζ2,c*ζ3
    return ζ0,ζ1,ζ2,ζ3 
end

#return π/6*@f(ρ_S)*∑(@f(x_S,k)*@f(d,k)^m for k ∈ @groups)

function ρ_S(model::SAFTgammaMieModel, V, T, z)
    return ρ_S(model.vrmodel,V,T,z)
end

function ζ_X_σ3(model::SAFTgammaMieModel, V, T, z,_d = @f(d_gc))
    vrmodel = model.vrmodel
    m = vrmodel.params.segment.values
    m̄ = dot(z, m)
    m̄inv = 1/m̄
    σ = model.params.sigma.values
    _mi = model.params.mixedsegment.values
    ρS = N_A/V*m̄
    _ζ_X = zero(first(_d))
    kρS = ρS* π/6
    σ3_x = _ζ_X
    r1 =_ζ_X
    for i ∈ @groups
        mi = _mi[i]
        x_Si = dot(mi,z)*m̄inv
        σ3_x += x_Si*x_Si*(σ[i,i]^3)
        _ζ_X += x_Si*x_Si*(_d[i])^3   
        for j ∈ 1:i-1
            mj = _mi[j]
            x_Sj = dot(mj,z)*m̄inv
            σ3_x += 2*x_Si*x_Sj*(σ[i,j]^3)
            dij = 0.5*(_d[i] + _d[j])
            r1 = 2*x_Si*x_Sj*dij^3          
            _ζ_X += r1
        end
    end

    #return π/6*@f(ρ_S)*∑(@f(x_S,i)*@f(x_S,j)*(@f(d,i)+@f(d,j))^3/8 for i ∈ comps for j ∈ comps)
    return kρS*_ζ_X,σ3_x
end

function σ3x(model::SAFTgammaMieModel, V, T, z)
    vrmodel = model.vrmodel
    m = vrmodel.params.segment.values
    m̄ = dot(z, m)
    m̄inv = 1/m̄
    σ = model.params.sigma.values
    _mi = model.params.mixedsegment.values
    σ3_x =  zero(first(z))
    for i ∈ @groups
        mi = _mi[i]
        x_Si = dot(mi,z)*m̄inv
        σ3_x += x_Si*x_Si*(σ[i,i]^3)
        for j ∈ 1:i-1
            mj = _mi[j]
            x_Sj = dot(mj,z)*m̄inv
            σ3_x += 2*x_Si*x_Sj*(σ[i,j]^3)
        end
    end
    return σ3_x
end

function  a_hs(model::SAFTgammaMieModel, V, T, z,_data = @f(data))
    _,vrdata = _data
    return a_hs(model.vrmodel,V,T,z,vrdata)
end

function a_chain(model::SAFTgammaMieModel, V, T, z,_data = @f(data))
    _,vrdata = _data
    return a_chain(model.vrmodel,V,T,z,vrdata)
end

function a_assoc(model::SAFTgammaMieModel, V, T, z,_data = @f(data))
    _,vrdata = _data
    return a_assoc(model.vrmodel,V,T,z,vrdata)
end

function a_disp(model::SAFTgammaMieModel, V, T, z,_data = @f(data))
    groups = @groups
    gcmodel = model
    model = gcmodel.vrmodel
    _d,vrdata = _data
    _,ρS,ζi,_ζ_X,_ζst,_ = vrdata
    comps = @comps
    l = length(comps)
    ∑z = ∑(z)
    m = model.params.segment.values
    _mi = gcmodel.params.mixedsegment.values
    _ϵ = gcmodel.params.epsilon.values
    _λr = gcmodel.params.lambda_r.values
    _λa = gcmodel.params.lambda_a.values
    _σ = gcmodel.params.sigma.values
    m̄ = dot(z, m)
    m̄inv = 1/m̄
    a₁ = zero(V+T+first(z))
    a₂ = a₁
    a₃ = a₁
    _ζst5 = _ζst^5
    _ζst8 = _ζst^8
    _KHS = @f(KHS,_ζ_X,ρS)
    for i ∈ groups
        j = i
        mi = _mi[i]
        x_Si = dot(mi,z)*m̄inv
        ϵ = _ϵ[i,j]
        λa = _λa[i,j]
        λr = _λr[i,j] 
        σ = _σ[i,j]
        _C = @f(Cλ,λa,λr)
        dij = _d[i]
        x_0ij = σ/dij
        dij3 = dij^3
        x_0ij = σ/dij
        #calculations for a1 - diagonal
        aS_1_a = @f(aS_1,λa,_ζ_X)
        aS_1_r = @f(aS_1,λr,_ζ_X)
        B_a = @f(B,λa,x_0ij,_ζ_X)
        B_r = @f(B,λr,x_0ij,_ζ_X)
        a1_ij = (2*π*ϵ*dij3)*_C*ρS*
        (x_0ij^λa*(aS_1_a+B_a) - x_0ij^λr*(aS_1_r+B_r))

        #calculations for a2 - diagonal
        aS_1_2a = @f(aS_1,2*λa,_ζ_X)
        aS_1_2r = @f(aS_1,2*λr,_ζ_X)
        aS_1_ar = @f(aS_1,λa+λr,_ζ_X)
        B_2a = @f(B,2*λa,x_0ij,_ζ_X)
        B_2r = @f(B,2*λr,x_0ij,_ζ_X)
        B_ar = @f(B,λr+λa,x_0ij,_ζ_X)
        α = _C*(1/(λa-3)-1/(λr-3))
        f1,f2,f3,f4,f5,f6 = @f(f123456,α)
        _χ = f1*_ζst+f2*_ζst5+f3*_ζst8
        a2_ij = π*_KHS*(1+_χ)*ρS*ϵ^2*dij3*_C^2 *
        (x_0ij^(2*λa)*(aS_1_2a+B_2a)
        - 2*x_0ij^(λa+λr)*(aS_1_ar+B_ar)
        + x_0ij^(2*λr)*(aS_1_2r+B_2r))
        
        #calculations for a3 - diagonal
        a3_ij = -ϵ^3*f4*_ζst * exp(f5*_ζst+f6*_ζst^2)
        #adding - diagonal
        a₁ += a1_ij*x_Si*x_Si
        a₂ += a2_ij*x_Si*x_Si
        a₃ += a3_ij*x_Si*x_Si
        for j ∈ 1:i-1
            mj = _mi[j]
            x_Sj = dot(mj,z)*m̄inv
            ϵ = _ϵ[i,j]
            λa = _λa[i,j]
            λr = _λr[i,j] 
            σ = _σ[i,j]
            _C = @f(Cλ,λa,λr)
            dij = 0.5*(_d[i]+_d[j])
            x_0ij = σ/dij
            dij3 = dij^3
            x_0ij = σ/dij
            #calculations for a1
            a1_ij = (2*π*ϵ*dij3)*_C*ρS*
            (x_0ij^λa*(@f(aS_1,λa,_ζ_X)+@f(B,λa,x_0ij,_ζ_X)) - x_0ij^λr*(@f(aS_1,λr,_ζ_X)+@f(B,λr,x_0ij,_ζ_X)))

            #calculations for a2
            α = _C*(1/(λa-3)-1/(λr-3))
            f1,f2,f3,f4,f5,f6 = @f(f123456,α)
            _χ = f1*_ζst+f2*_ζst5+f3*_ζst8
            a2_ij = π*_KHS*(1+_χ)*ρS*ϵ^2*dij3*_C^2 *
            (x_0ij^(2*λa)*(@f(aS_1,2*λa,_ζ_X)+@f(B,2*λa,x_0ij,_ζ_X))
            - 2*x_0ij^(λa+λr)*(@f(aS_1,λa+λr,_ζ_X)+@f(B,λa+λr,x_0ij,_ζ_X))
            + x_0ij^(2*λr)*(@f(aS_1,2λr,_ζ_X)+@f(B,2*λr,x_0ij,_ζ_X)))
            
            #calculations for a3
            a3_ij = -ϵ^3*f4*_ζst * exp(f5*_ζst+f6*_ζst^2)
            #adding
            a₁ += 2*a1_ij*x_Si*x_Sj
            a₂ += 2*a2_ij*x_Si*x_Sj
            a₃ += 2*a3_ij*x_Si*x_Sj            
        end
    end
    a₁ = a₁*m̄/T/∑z
    a₂ = a₂*m̄/(T*T)/∑z
    a₃ = a₃*m̄/(T*T*T)/∑z
    #@show (a₁,a₂,a₃)
    adisp =  a₁ + a₂ + a₃ 
    return adisp
end

const SAFTγMieconsts =(
    A = SA[ 0.81096    1.7888   -37.578   92.284;
          1.02050  -19.341    151.26  -463.50 ;
         -1.90570   22.845   -228.14   973.92 ;
          1.08850   -6.1962   106.98  -677.64  ],

    ϕ = [[7.5365557, -359.440,  1550.9, -1.199320, -1911.2800,  9236.9  ],
        [-37.604630,  1825.60, -5070.1,  9.063632,  21390.175, -129430 ],
        [71.745953, -3168.00,  6534.6, -17.94820, -51320.700,  357230  ],
        [-46.835520,  1884.20, -3288.7,  11.34027,  37064.540, -315530 ],
        [-2.4679820,- 0.82376, -2.7171,  20.52142,  1103.7420,  1390.2 ],
        [-0.5027200, -3.19350,  2.0883, -56.63770, -3264.6100, -4518.2 ],
        [8.0956883,  3.70900,  0.0000,  40.53683,  2556.1810,  4241.6  ]],

    u = [0.26356031971814109102031,1.41340305910651679221800,3.59642577104072208122300,7.08581000585883755692200,12.6408008442757826594300],
    w = [0.5217556105828086524759,0.3986668110831759274500,7.5942449681707595390e-2,3.6117586799220484545e-3,2.3369972385776227891e-5],

         c  = [0.0756425183020431	-0.128667137050961	 0.128350632316055	-0.0725321780970292	   0.0257782547511452  -0.00601170055221687	  0.000933363147191978  -9.55607377143667e-05  6.19576039900837e-06 -2.30466608213628e-07 3.74605718435540e-09
                0.134228218276565	    -0.182682168504886 	 0.0771662412959262	-0.000717458641164565 -0.00872427344283170	0.00297971836051287	 -0.000484863997651451	 4.35262491516424e-05 -2.07789181640066e-06	4.13749349344802e-08 0
               -0.565116428942893	     1.00930692226792   -0.660166945915607	 0.214492212294301	  -0.0388462990166792	0.00406016982985030	 -0.000239515566373142	 7.25488368831468e-06 -8.58904640281928e-08	0	                 0
               -0.387336382687019	    -0.211614570109503	 0.450442894490509	-0.176931752538907	   0.0317171522104923  -0.00291368915845693	  0.000130193710011706  -2.14505500786531e-06  0	                0	                 0
                2.13713180911797	    -2.02798460133021 	 0.336709255682693	 0.00118106507393722  -0.00600058423301506	0.000626343952584415 -2.03636395699819e-05	 0	                   0	                0	                 0
               -0.300527494795524	     2.89920714512243   -0.567134839686498	 0.0518085125423494	  -0.00239326776760414	4.15107362643844e-05  0	                     0	                   0	                0                    0
               -6.21028065719194	    -1.92883360342573	 0.284109761066570	-0.0157606767372364	   0.000368599073256615	0 	                  0	                     0	                   0	                0	                 0
                11.6083532818029	     0.742215544511197  -0.0823976531246117	 0.00186167650098254   0	                0	                  0	                     0	                   0	                0	                 0
               -10.2632535542427	    -0.125035689035085	 0.0114299144831867	 0	                   0	                0	                  0	                     0	                   0	                0	                 0
                4.65297446837297	    -0.00192518067137033 0	                 0	                   0	                0	                  0	                     0	                   0	                0	                 0
               -0.867296219639940	     0	                 0	                 0	                   0	                0	                  0	                     0	                   0	                0	                 0],
)
