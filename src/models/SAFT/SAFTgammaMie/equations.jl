function x0_volume_liquid(model::SAFTgammaMieModel,T,z)
    v_lb = lb_volume(model,z)
    return v_lb*2.0
end

function lb_volume(model::SAFTgammaMieModel, z)
    vk  = model.groups.n_flattenedgroups
    seg = model.params.segment.values
    S   = model.params.shapefactor.values
    σ = model.params.sigma.values
    lb_v = zero(Base.promote_eltype(model,z))
    for i in @comps
        vki = vk[i]
        lb_vi = zero(lb_v)
        for k in @groups(i)
            lb_vi += vki[k]*seg[k]*S[k]*σ[k,k]^3
        end
        lb_v += z[i]*lb_vi
    end

    return π/6*N_A*lb_v
    return val
end

function T_scale(model::SAFTgammaMieModel,z)
    return T_scale(model.vrmodel,z)
end

function T_scales(model::SAFTgammaMieModel)
    return T_scales(model.vrmodel)
end

function p_scale(model::SAFTgammaMieModel,z)
    V = zero(first(z))
    T = zero(first(z))
    σ̄3 = @f(σ3x)
    ϵ̄ = T_scale(model,z)
    val    = σ̄3*N_A/R̄/ϵ̄
    return 1/val
end

getsites(model::SAFTgammaMieModel) = model.vrmodel.sites
assoc_shape(model::SAFTgammaMieModel) = assoc_shape(model.vrmodel)
function a_res(model::SAFTgammaMieModel, V, T, z)
    _data = @f(data)
    dgc,X,vrdata = _data
    _,ρS,ζi,_ζ_X,_ζst,σ3x,m̄ = vrdata
    vrdata_disp = (dgc,ρS,ζi,_ζ_X,_ζst,σ3x,m̄)
    return @f(a_hs,_data) + a_disp(model,V,T,X,vrdata_disp)/sum(z) + @f(a_chain,_data) + @f(a_assoc,_data)
end

function a_mono(model::SAFTgammaMieModel, V, T, z,_data = @f(data))
    dgc,X,vrdata = _data
    _,ρS,ζi,_ζ_X,_ζst,σ3x,m̄ = vrdata
    vrdata_disp = (dgc,ρS,ζi,_ζ_X,_ζst,σ3x,m̄)
    return @f(a_hs,_data) + a_disp(model,V,T,X,vrdata_disp)/sum(z)
end

function data(model::SAFTgammaMieModel, V, T, z)
    m̄ = dot(z, model.vrmodel.params.segment.values)
    X = @f(X_gc)
    _d_gc = d(model,V,T,X)
    _d_gc_av = @f(d_gc_av,_d_gc)
    ζi = ζ0123(model,V,T,X,_d_gc)
    _ζ_X,σ3x = ζ_X_σ3(model,V,T,X,_d_gc,m̄)
    _ρ_S = N_A/V*m̄
    _ζst = _ζst = σ3x*_ρ_S*π/6
    vrdata = (_d_gc_av,_ρ_S,ζi,_ζ_X,_ζst,σ3x,m̄)
    return (_d_gc,X,vrdata)
end

function packing_fraction(model::SAFTgammaMieModel,_data::Tuple)
    _,_,vrdata = _data
    return packing_fraction(model.vrmodel,vrdata)
end

function X_gc(model::SAFTgammaMieModel,V,T,z)
    mi  = group_matrix(model.groups)::Matrix{Float64}
    mm = model.params.segment.values
    TT = Base.promote_eltype(mm,1.0,z)
    X = Vector{TT}(undef,length(model.groups.flattenedgroups))
    mul!(X,mi,z)
    X ./= mm
    return X
end

function d_gc_av(model::SAFTgammaMieModel,V,T,z,_d_gc = d(model,V,T,@f(X_gc)))
    _d = zeros(eltype(_d_gc),length(z))
    _0 = zero(eltype(_d))
    _zz = model.groups.n_groups_cache
    for i ∈ @comps
        _z = _zz[i]
        ∑zinv2 = 1/(sum(_z)^2)
        di = _0
        for k ∈ @groups
            zk = _z[k]
            iszero(zk) && continue
            dk = _d_gc[k]
            di += zk*zk*dk^3
            for l ∈ 1:k - 1
                di += 0.25*zk*_z[l]*(dk + _d_gc[l])^3 # 2*(di + dj/2)^3 = di 
            end
        end
        _d[i] = cbrt(di*∑zinv2)
    end
    #d̂[i] = cbrt(∑(∑(@f(ẑ,i,k)*@f(ẑ,i,l)*@f(d,k,l)^3 for l ∈ @groups) for k ∈ @groups))
    return _d
end

#return π/6*@f(ρ_S)*∑(@f(x_S,k)*@f(d,k)^m for k ∈ @groups)

function σ3x(model::SAFTgammaMieModel, V, T, z)
    X = @f(X_gc)
    m = model.params.segment.values
    m̄ = dot(z, model.vrmodel.params.segment.values)
    m̄inv = 1/m̄
    σ = model.params.sigma.values
    σ3_x = zero(first(z))
    for i ∈ @groups
        x_Si = X[i]*m[i]*m̄inv
        σ3_x += x_Si*x_Si*(σ[i,i]^3)
        for j ∈ 1:i-1
            x_Sj = X[j]*m[j]*m̄inv
            σ3_x += 2*x_Si*x_Sj*(σ[i,j]^3)
        end
    end
    return σ3_x
end

function  a_hs(model::SAFTgammaMieModel, V, T, z,_data = @f(data))
    _,_,vrdata = _data
    return a_hs(model.vrmodel,V,T,z,vrdata)
end

function a_chain(model::SAFTgammaMieModel, V, T, z,_data = @f(data))
    _,_,vrdata = _data
    return a_chain(model.vrmodel,V,T,z,vrdata)
end


function a_assoc(model::SAFTgammaMieModel, V, T, z,_data = @f(data))
    _,_,vrdata = _data
    return a_assoc(model.vrmodel,V,T,z,vrdata)
end

function Δ(model::SAFTgammaMieModel, V, T, z, i, j, a, b,_data = @f(data))
    vrdata = _data[3]
    return Δ(model.vrmodel,V,T,z,i,j,a,b,vrdata)
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

########
#=
Optimizations for single component SAFTgammaMieModel
=#

#######

function d_gc_av(model::SAFTgammaMieModel,V,T,z::SingleComp,_d_gc = d(model,V,T,@f(X_gc)))
    _0 = zero(eltype(_d_gc))
    _z = only(model.groups.n_groups_cache)  
    ∑zinv2 = 1/(sum(_z)^2)
    di = _0
    for k ∈ @groups
        zk = _z[k]
        iszero(zk) && continue
        dk = _d_gc[k]
        di += zk*zk*dk^3
        for l ∈ 1:k - 1
            di += 0.25*zk*_z[l]*(dk + _d_gc[l])^3 # 2*(di + dj/2)^3 = di 
        end
    end
    return SA[cbrt(di*∑zinv2)]
end

#=
function X_homotopy(model,V,T,z,_data = Clapeyron.data(model,V,T,z))
    _Δ = Clapeyron.__delta_assoc(model,V,T,z,_data)
    n = Clapeyron.assoc_pair_length(model)
    K = Clapeyron.assoc_site_matrix(model,V,T,z,_data,_Δ)
    sitesparam = Clapeyron.getsites(model)
    idxs = sitesparam.n_sites.p
    nn = length(sitesparam.n_sites.v)
    @var xi[1:nn]
    system = System(K*xi .* xi + xi - ones(Int,nn))
    result = solve(system) 
    real_result = real_solutions(result)
    function f(xx)
        all(xk -> 0 <= xk <= 1, xx)
    end
    display(real_result)
    idx = findfirst(f,real_result)
    xsol = real_result[idx]
    
    return Clapeyron.PackedVofV(idxs,xsol)
end
=#