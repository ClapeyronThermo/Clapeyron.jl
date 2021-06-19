#aproximates liquid volume at a known pressure and t,
#by using isothermal compressibility
function volume_compress(model,p,T,z=SA[1.0];V0=x0_volume(model,p,T,z,phase=:liquid),max_iters=100)
    logV0 = log(V0)
    function f_fixpoint(_V)
        _V = exp(_V)
        _p,dpdV = p∂p∂V(model,_V,T,z)
        β = -1/_V*dpdV^-1
        _Δ =  -(p-_p)*β
        sign_Δ = sign(_Δ)
        Δ = abs(_Δ)
        Vv = _V*exp(sign_Δ*Δ^(1-_Δ))#^((1+Δ)^4)
        return log(Vv)
    end
        res = Solvers.fixpoint(f_fixpoint,logV0,Solvers.SimpleFixPoint(),rtol = 1e-12,max_iters=max_iters)
        return exp(res)
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
        #return value at spinodal
        return -2*B
    end
    #only the left root has physical meaning
    return (-b + sqrt(b*b-4*a*c))/(2*a)
end

#(z = pV/RT)
#(RT/p = V/z)

function volume(model::EoSModel,p,T,z=SA[1.0];phase=:unknown,threaded=true)

    fp(_V) = log(pressure(model,_V,T,z)/p)

    el = false
    ev = false
    Vg = NaN
    Vl = NaN
#Threaded version
    phase = Symbol(phase)
    if phase != :unknown
        V0 = x0_volume(model,p,T,z,phase=phase)
        #return Solvers.ad_newton(fp,Vg0)
        return volume_compress(model,p,T,z,V0=V0)
    end

    if threaded
        try
            Vg0 = x0_volume(model,p,T,z,phase=:v)
            Vg = fetch(Threads.@spawn volume_compress(model,$p,$T,$z;V0=Vg0))
        catch 
            ev  = true
        end

        try
            Vl0 = x0_volume(model,p,T,z,phase=:l)
            Vl = fetch(Threads.@spawn volume_compress(model,$p,$T,$z;V0=$Vl0))
        catch 
            el  = true
        end
    else
        Vg0 = x0_volume(model,p,T,z,phase=:v)
        Vl0 = x0_volume(model,p,T,z,phase=:l)
        try
            Vg =  volume_compress(model,p,T,z,V0=Vg0)
        catch 
            ev=true
        end
        #Vg = Solvers.ad_newton(fp,Vg0)
        try
            Vl =  volume_compress(model,p,T,z,V0=Vl0)
        catch 
            el=true
        end
    end

# Serial version
#=
    Vg0 = volume_virial(model,p,T,z)
    Vg =Solvers.ad_newton(fp,Vg0,rtol=1e-08)

    Vl =  volume_compress(model,p,T,z)
=#
    function gibbs(V)
        _df,f =  ∂f(model,V,T,z)
        dV,dt = _df
        if abs((p+dV)/p) > 0.03
            return Inf
        else
            return f  +p*V
        end
    end
    #this catches the supercritical phase as well

    if isnan(Vl) || el
        return Vg
    end
    if isnan(Vg) || ev
        return Vl
    end
    if Vl ≈ Vg
        return Vl
    end

    if ev && el
        error("Failed to converge to a root")
    end
        gg = gibbs(Vg)
        gl = gibbs(Vl)
        #@show Vg,Vl
        #@show gg,gl
        return ifelse(gg<gl,Vg,Vl)

end

export volume
