#aproximates liquid volume at a known pressure and t,
#by using isothermal compressibility
function volume_compress(model,p,T,z=SA[1.0])
    v0 = vcompress_v0(model,p,T,z)
    function f_fixpoint(_v)
        _p,dpdv = p∂p∂v(model,_v,T,z)
        β = -1/_v*dpdv^-1
        _Δ =  -(p-_p)*β
        sign_Δ = sign(_Δ)
        Δ = abs(_Δ)
        vv = _v*exp(sign_Δ*Δ^(1-_Δ))#^((1+Δ)^4)
        return vv
    end
    return Solvers.fixpoint(f_fixpoint,v0,Solvers.SimpleFixPoint(),rtol = 1e-12)
    #return Roots.find_zero(f_fixpoint,v0)
end

function vcompress_v0(model,p,T,z=SA[1.0])
    lb_v   = lb_volume(model,z)
    v0 = 1.1*lb_v
    return v0
end

function volume_virial(model,p,T, z=SA[1.] )
    B = second_virial_coefficient(model,T,z)
    a = p/(R̄*T)
    b = -1
    c = -B
    Δ = b*b-4*a*c
    n = sum(z)
    if Δ <= 0
        #virial approximation could not be calculated
        #degrade to ideal aprox
        return n*R̄*T/p 
        
    end 
    return (-b + sqrt(b*b-4*a*c))/(2*a)
end

#(z = pv/rt)
#(RT/p = v/z)

function volume(model::EoSModel,p,T,z=SA[1.0];phase="unknown")

    fp(_v) = pressure(model,_v,T,z) - p

#Threaded version

    if is_liquid(phase)
        return volume_compress(model,p,T,z)
    elseif is_vapour(phase) | is_supercritical(phase)
        vg0 = volume_virial(model,p,T,z)
        return Solvers.ad_newton(fp,vg0,rtol=1e-08)
    else #unknown handled here
    end
    
    _vg = Threads.@spawn begin
        vg0 = volume_virial(model,$p,$T,$z)
        Solvers.ad_newton(fp,vg0,rtol=1e-08)
    end
        
    _vl = Threads.@spawn volume_compress(model,$p,$T,$z)
    #fp(_v) = pressure(model,_v,T,z) - p
    vg = fetch(_vg)
    vl = fetch(_vl)


# Serial version
#=
    vg0 = volume_virial(model,p,T,z)
    vg =Solvers.ad_newton(fp,vg0,rtol=1e-08)

    vl =  volume_compress(model,p,T,z)
=#
    function gibbs(v)
        _df,f =  ∂f(model,v,T,z)
        dv,dt = _df
        if abs((p+dv)/p) > 0.03
            return Inf
        else
            return f  +p*v
        end
    end
    #this catches the supercritical phase as well
    if vl ≈ vg
        return vl
    end  
        gg = gibbs(vg)
        gl = gibbs(vl)
        #@show vg,vl
        #@show gg,gl
        return ifelse(gg<gl,vg,vl)
    
end
