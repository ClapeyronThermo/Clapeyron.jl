# ## Old pressure solver


# function volume(model::EoSModel, p, T,  z=SA[1.]; phase = "unknown")
#     N = length(p)

#     ub = [Inf]
#     lb = lb_volume(model,z; phase = phase)
    
#     x0 = x0_volume(model,p,T,z; phase = phase)
#     f = v -> eos(model, exp10(v[1]), T,z) + exp10(v[1])*p
#     #looking at best phase using tunneling
#     if phase == "unknown"
#         (f_best,v_best) = Solvers.tunneling(f,lb,ub,x0)
#         #@show eos.(Ref(model),exp10.(v_best),T,Ref(z)) .- p.*exp10.(v_best)

#         return exp10(v_best[1])
#     else
#         _ub = 100*T/p
#         _lb = exp10(only(lb))
#         f0(vx) = pressure(model,vx,T,z,phase=phase) - p        
#         _v0 = exp10(only(x0))
#         #try direct newton solving
#         vobj = Solvers.ad_newton(f0,_v0)
#         if (_lb <= vobj <= _ub)
#             return vobj
#         else
#             vobj = Roots.find_zero(f0,(_lb,_ub),FalsePosition())
#         end
#         #=
#         opt_min = NLopt.Opt(:LD_MMA, length(ub))
#         opt_min.lower_bounds = lb
#         opt_min.upper_bounds = ub
#         opt_min.xtol_rel     = 1e-8
#         obj_f0 = x -> f(x)
#         obj_f  = (x,g) -> Solvers.NLopt_obj(obj_f0,x,g)
#         opt_min.min_objective =  obj_f
#         (f_min,v_min) = NLopt.optimize(opt_min, x0)
#         #@show eos.(Ref(model),exp10.(v_min),T,Ref(z)) .- p.*exp10.(v_min)
#         return exp10(v_min[1])
#         =#
#     end
# end

# function Jac_Sat(model::EoSModel, J, T, v_l, v_v)
#     #components = model.components
#     fun(x) = eos(model, x[2], T,SA[x[1]])
#     d2f(x) = ForwardDiff.hessian(fun,x)
#     d2f_l = d2f(SA[one(v_l),v_l])
#     d2f_v = d2f(SA[one(v_v),v_v])
#     (p_scale,μ_scale) = scale_sat_pure(model)
#     J[1,1] =  v_l[1]*d2f_l[2,2]*log(10)*p_scale
#     J[1,2] = -v_v[1]*d2f_v[2,2]*log(10)*p_scale
#     J[2,1] =  v_l[1]*d2f_l[1,2]*log(10)*μ_scale
#     J[2,2] = -v_v[1]*d2f_v[1,2]*log(10)*μ_scale
# end

# function Obj_Jac_sat(model::EoSModel, F, J, T, v_l, v_v)
#     Obj_Sat(model, F, T, v_l, v_v)
#     Jac_Sat(model, J, T, v_l, v_v)
#     F, J
# end

# function Jvop_sat(x,model::EoSModel,T)
#     function Jac_satV(Fv, v)
#         #components = model.components
#         fun(x) = eos(model, x[2], T,[x[1]])
#         d2f(z) = ForwardDiff.hessian(fun,z)
#         d2f_l = d2f([1,x[1]])
#         d2f_v = d2f([1,x[2]])
#         (p_scale,μ_scale) = scale_sat_pure(model)
#         Fv[1,] =  x[1]*d2f_l[2,2]*log(10)*p_scale*v[1]-x[2]*d2f_v[2,2]*log(10)*p_scale*v[2]
#         Fv[2,] =  x[1]*d2f_l[1,2]*log(10)*μ_scale*v[1]-x[2]*d2f_v[1,2]*log(10)*μ_scale*v[2]
#     end
#     LinearMap(Jac_satV, length(x))
# end




# ## Mixture critical point solver
# # function crit_mix(model::SAFT,x_c)
# #     components = model.components
# #     z  = create_z(model,x_c)
# #     f! = (F,x) -> Obj_Crit_mix(model, F, exp10(x[2]), x[1]*prod(model.params.epsilon[i]^z[i] for i in components), x_c)
# #     x0 = [1.5,log10(π/6*N_A*sum(z[i]*model.params.segment[i]*model.params.sigma[i]^3 for i in components)/0.15)]
# #     r  = Solvers.nlsolve(f!,x0)
# #     T_c = r.zero[1]*prod(model.params.epsilon[i]^z[i] for i in components)
# #     v_c = exp10(r.zero[2])
# #     p_c = pressure(model, v_c, T_c, x_c)
# #     return (T_c, p_c, v_c)
# # end
# # #
# # function Obj_Crit_mix(model::SAFT, F, v_c,T_c,x_c)
# #     fun(x)  = eos(model, create_z(model, [x[1],1-x[1]]), v_c, T_c)
# #     df(x)   = ForwardDiff.derivative(fun,x)
# #     d2f(x)  = ForwardDiff.derivative(df,x)
# #     d3f(x)  = ForwardDiff.derivative(d2f,x)
# #     F[1] = d2f(x_c[1])
# #     F[2] = d3f(x_c[1])
# #     println(F)
# #     println(v_c)
# #     println(T_c)
# # end

# function Jac_bubble_pressure(model, J, T, v_l, v_v, x, y)
#     components = model.components
#     append!(y,1-sum(y[i] for i in 1:(length(components)-1)))

#     fun(z) = eos(model, z[1:end-1], z[end], T)
#     df(z)  = ForwardDiff.gradient(fun,z)
#     d2f(z) = ForwardDiff.hessian(fun,z)
#     X = deepcopy(x)
#     Y = deepcopy(y)
#     d2f_l = d2f(append!(X,v_l))
#     d2f_v = d2f(append!(Y,v_v))
#     for i in 1:(length(components))
#         J[i,1] =  log(10)*v_l*d2f_l[i,end]/R̄/model.params.epsilon[components[i]]
#         J[i,2] = -log(10)*v_v*d2f_v[i,end]/R̄/model.params.epsilon[components[i]]
#     end
#     J[end,1] = log(10)*v_l*d2f_l[end,end]*model.params.sigma[components[1]]^3*N_A/R̄/model.params.epsilon[components[1]]
#     J[end,2] =-log(10)*v_v*d2f_v[end,end]*model.params.sigma[components[1]]^3*N_A/R̄/model.params.epsilon[components[1]]

#     for j in 1:(length(components)-1)
#         J[end,j+2] = (d2f_v[end,end-1]-d2f_v[end,j])*model.params.sigma[components[1]]^3*N_A/R̄/model.params.epsilon[components[1]]
#     end

#     for i in 1:(length(components))
#         for j in 1:(length(components)-1)
#                 J[i,j+2]= -(d2f_v[i,j]-d2f_v[i,end-1])/R̄/model.params.epsilon[components[i]]
#         end
#     end
# end

# ##Old sat pure code in case of points near the critical point
# # if debug
#         #     if error_val[] !== nothing
#         #         @warn "initial saturation calculation failed with error $error_val[]"
#         #     end
#         # end
#         # result[] = res0
#         # #convergence not achieved, trying critical aproximation


#         # T_c,P_c,V_c = crit_pure(model)
#         # ΔT = (T_c - T[i])
#         # ΔT <= 8*eps(ΔT) && throw(DomainError(T[i],"input temperature $T is too close or higher than critical temperature of the model $T_c"))
#         # Tr  = T[i]/T_c
#         # V0 = x0_sat_pure_crit(model,T[i],T_c,P_c,V_c)

#         # try_sat_pure(model,V0,f!,T[i],result,error_val,converged)
#         # if converged[]
#         #     append!(p_sat,result[][1])
#         #     append!(V_l,result[][2])
#         #     append!(V_v,result[][3])
#         #     V0 = log10.([V_l[i],V_v[i]])
#         # else
#         #     @warn "the procedure converged to a trivial value at T=$T"
#         #     return result[]
#         # end


#         # if debug
#         #     throw(error_val[])
#         # else
#         #     throw("unable to calculate equilibria at T=$T")
#         # end
    