function rr_bubble_pressure_refine(model,x,y,vl,vv,T)
    pl = pressure(model,vl,T,x)
    pv = pressure(model,vv,T,y)
    p = (pl+pv)/2
    vl = volume(model,p,T,x,phase=:l)
    vv = volume(model,p,T,y,phase=:v)
    μ_l = vt_chemical_potential(model,v_l,T,x)
    #@show "c"
    μ_v = vt_chemical_potential(model,v_v,T,y)
    K = μ_v/μ_l
    y = K.*x
    @show vl,vv,y
    return model,x,y,vl,vv,T
end

    
