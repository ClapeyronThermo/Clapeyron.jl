#=
SPUNG = State Research Program for Utilization of Natural Gas

based on:
Unification of the two-parameter equation of state and the principle ofcorresponding states
Jørgen Mollerup

the idea is to use mollerup shape factors and PropaneRef as the EoS
but any function that provides a_scaled and shape_factors will suffice
=#

struct SPUNG{S<:EoSModel,E<:EoSModel} <: EoSModel
    shape_model::S
    shape_ref::S
    model_ref::E
end

function SPUNG(components::Vector{String},refmodel=PropaneRef(),shapemodel::SHAPE=SRK(components)) where SHAPE<:EoSModel
    refname = component_names(refmodel)
    shape_ref = SHAPE.name.wrapper(refname)
    return SPUNG(shapemodel,shape_ref,refmodel)
end


function eos(model::SPUNG,V,T,z=SA[1.0],phase="unknown")
    f,h = shape_factors(model,V,T,z)
    T0 = T/f
    V0 = V/h
    return eos(model.model_ref,V0,T0)
end

function shape_factors(model::SPUNG{<:ABCubicModel},V,T,z=SA[1.0])
    n = sum(z)
    x = z * (1/n)
    a,b = cubic_ab(model.shape_model,T,x)
    a0,b0 = cubic_ab(model.shape_ref,T,SA[1.0])
    h = b/b0
    fh = n*a/a0
    f = fh/h
    return f,h
end

mw(model::SPUNG) = mw(model.shape_model)

function Base.show(io::IO,mime::MIME"text/plain",model::SPUNG)
    println(io,"Extended Corresponding States model")
    println(io," reference model: ",model.model_ref)
    print(io," shape model: ",model.shape_model)
end

function Base.show(io::IO,model::SPUNG)
    print(io,"SPUNG(",model.shape_ref,",",model.model_ref,")")
end

function lb_volume(model::SPUNG,z=SA[1.0];phase=:l)
    lb_v0 = lb_volume(model.model_ref)
    T0 = T_scale(model.model_ref)
    f,h = shape_factors(model,lb_v0,T0,z) #h normaly should be independent of temperature
    return lb_v0*h
end

function T_scale(model::SPUNG,z=SA[1.0])
    lb_v0 = lb_volume(model.model_ref)
    T0 = T_scale(model.model_ref)
    f,h = shape_factors(model,lb_v0,T0,z) #h normaly should be independent of temperature
    return T0*f
end

function p_scale(model::SPUNG,z=SA[1.0])
    lb_v0 = lb_volume(model.model_ref)
    T0 = T_scale(model.model_ref)
    p0 = p_scale(model.model_ref)
    f,h = shape_factors(model,lb_v0,T0,z) #h normaly should be independent of temperature
    ps = p0*f/h 
    return ps
end

#=
ideally we could perform SPUNG only providing x0, but i cant find the error here
function x0_sat_pure(model::SPUNG,T,z=SA[1.0])
    lb_v = lb_volume(model,z)
    f,h = shape_factors(model,lb_v,T,z)
    T0 = T/f 
    @show T0
    vl0,vv0 = exp10.(x0_sat_pure(model.model_ref,T0,SA[1.0]))
    @show vl = vl0*h
    @show vv = vv0*h
    return [log10(vl),log10(vv)]
end
=#

#overloading sat_pure for SPUNG directly seems to be the way
function sat_pure(model::SPUNG,T::Real)
    lb_v = lb_volume(model,SA[1.0])
    f,h = shape_factors(model,lb_v,T,SA[1.0])
    T0 = T/f
    psat0,vl0,vv0 = sat_pure(model.model_ref,T0)
    p = pressure(model,vv0*h,T)
     return (p,vl0*h,vv0*h)
end

function sat_pure_p(model::SPUNG,p::Real)
    return naive_sat_pure_p(model,p)
end


function general_shape_factors(model::SPUNG,V,T,z=SA[1.0])
    n = sum(z)
    x = z * (1/n)
    RT = R̄*T
    b = lb_volume(model.shape_model,x)
    b0 = lb_volume(model.shape_ref,x)
    B = second_virial_coefficient(model.shape_model,T,x)
    B0 = second_virial_coefficient(model.shape_ref,T)
    #B = b-a/RT
    #a/RT = b-B
    #a = RT(b-B)
    a = RT*(b-B)
    a0 = RT*(b0-B0)
    tau = 1/(1-4(B/b)) 
    tau0 = 1(1-4(B0/b0))
    Tc = T_scale(model.shape_model,x)
    Tc0 = T_scale(model.shape_ref,x)

    f0 = (tau*Tc)/(tau0*Tc0)
    #f0 = tau/tau0
    #@show T/f0
     #T0 = Roots.find_zero(f0_f,T/f0)
    #@show T0
     #B0 = second_virial_coefficient(model.shape_ref,T0)
    #tau0 = 1-(B0/b0) 
    #tau0 = 1-(B0/b0)
    #f = tau/tau0
    
    #a,b = cubic_ab(model.shape_model,T,x)
    #a0,b0 = cubic_ab(model.shape_ref,T,SA[1.0])
    h = b/b0
    #fh = n*a/a0
    #f = fh/h
    f = f0
    return f,h
end

#=
Tc = 0.26+2.1R
R = λ-1


=#

export SPUNG