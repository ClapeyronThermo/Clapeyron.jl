"""
SPUNG: State Research Program for Utilization of Natural Gas

based on:
Unification of the two-parameter equation of state and the principle ofcorresponding states
Jørgen Mollerup

The idea is to use a "shape model" that provides a corresponding states parameters
and a "reference model" that implements a helmholtz energy function.

this SPUNG by default uses the propane reference equation of state (`PropaneRef`) as the reference model
and SRK for the shape model.
"""
struct SPUNG{S<:EoSModel,E<:EoSModel} <: EoSModel
    shape_model::S
    shape_ref::S
    model_ref::E
end

function SPUNG(components::Vector{String},refmodel=PropaneRef(),shapemodel=SRK(components),shaperef = SRK(component_names(refmodel)))
    model = SPUNG(shapemodel,shaperef,refmodel)
    return model
end


function eos(model::SPUNG,V,T,z=SA[1.0],phase=:unknown)
    f,h = shape_factors(model,V,T,z)
    T0 = T/f
    V0 = V/h
    return eos(model.model_ref,V0,T0)
end

function eos_res(model::SPUNG,V,T,z=SA[1.0],phase=:unknown)
    f,h = shape_factors(model,V,T,z)
    T0 = T/f
    V0 = V/h
    return eos_res(model.model_ref,V0,T0)
end

function shape_factors(model::SPUNG{<:ABCubicModel},V,T,z=SA[1.0])
    a,b = cubic_ab(model.shape_model,V,T,z)
    a0,b0 = cubic_ab(model.shape_ref,V,T,z)
    h = b/b0
    fh = a/a0
    f = fh/h
    return f,h
end

mw(model::SPUNG) = mw(model.shape_model)
molecular_weight(model::SPUNG,z=SA[1.0]) = molecular_weight(model.shape_model,z)

function Base.show(io::IO,mime::MIME"text/plain",model::SPUNG)
    println(io,"Extended Corresponding States model")
    println(io," reference model: ",model.model_ref)
    print(io," shape model: ",model.shape_model)
end

function Base.show(io::IO,model::SPUNG)
    print(io,"SPUNG(",model.shape_ref,",",model.model_ref,")")
end

function lb_volume(model::SPUNG,z=SA[1.0])
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

# function p_scale(model::SPUNG,z=SA[1.0])
#     lb_v0 = lb_volume(model.model_ref)
#     T0 = T_scale(model.model_ref)
#     p0 = p_scale(model.model_ref)
#     f,h = shape_factors(model,lb_v0,T0,z) #h normaly should be independent of temperature
#     ps = p0*f/h
#     return ps
# end

#=
ideally we could perform SPUNG only providing x0, but i cant find the error here

=#
#overloading saturation_pressure for SPUNG directly seems to be the way
function saturation_pressure(model::SPUNG,T::Real,v0=[zero(T)/zero(T),zero(T)/zero(T)])
    lb_v0 = lb_volume(model.model_ref)
    f,h = shape_factors(model,lb_v0,T,SA[1.0]) #h normaly should be independent of temperature
    T0 = T/f
    if isnan(v0[1]) && isnan(v0[2])
        v0 = x0_sat_pure(model.model_ref,T0)
    else
        vl = exp10(v0[1])
        vv = exp10(v0[2])
        v0 = [log10(vl*h),log10(vv*h)]
    end
    psat0,vl0,vv0 = saturation_pressure(model.model_ref,T0,v0)
    p = pressure(model,vl0*h,T)
    return (p,vl0*h,vv0*h)
end

#============
uncomment when saturation_pressure_p is ready
====================
function saturation_pressure_p(model::SPUNG,p::Real)
    return naive_saturation_pressure_p(model,p)
end
=#

# function split_model(model::SPUNG)
#     #only the shape model is splittable
#     shape_model = split_model(model.shape_model)
#     len = length(shape_model)
#     shape_ref = fill(model.shape_ref,len)
#     model_ref = fill(model.model_ref,len)
#     return SPUNG.(shape_model,shape_ref,model_ref)
# end

function shape_factors(model::SPUNG,V,T,z=SA[1.0])
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
    tau0 = 1/(1-4(B0/b0))
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
