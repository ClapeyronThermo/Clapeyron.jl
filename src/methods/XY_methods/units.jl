#Units handling

# General unit stripping
@inline ustrip(x,f::F) where F = ustrip(unit_system(x),x,f)
@inline ustrip(::Nothing,x,f::F) where F = x

#unit stripping for compositions: must handle total compositions (mass and molar). TODO: how to mark a vector to be with "mass fraction" units?
@inline uzstrip(model,z::Number) = uzstrip(unit_system(z),model,SA[z])
@inline uzstrip(model,z) = uzstrip(unit_system(z),model,z)
@inline uzstrip(::Nothing,model,z) = z

#unit stripping for volume: must handle volume (total, molar, mass) and density (molar, mass)
@inline uvstrip(model,v,z) = uvstrip(unit_system(v),model,v,z)
@inline uvstrip(::Nothing,model,v,z) = v

#entropy/internal energy: must handle total, molar and mass variants
@inline uhstrip(model,h,z) = uhstrip(unit_system(h),model,h,z)
@inline uhstrip(::Nothing,model,h,z) = h

#entropy: must handle total, molar and mass variants
@inline usstrip(model,s,z) = usstrip(unit_system(s),model,s,z)
@inline usstrip(::Nothing,model,s,z) = s

#any units package needs to overload:
function unit_type end #unit_type(us::MyUnitSystem,output)
function solve_unit end #solve_unit(unit_sys,output1,output2)
#with_output_unit(res,t::Tuple{X,Y}) where {X,Y}
unit_system(x::AbstractArray{T}) where T = unit_system(T)
unit_system(x) = nothing
unit_system(x::T,::Nothing) where T = x
unit_system(::Nothing,y::T) where T = y
unit_system(::Nothing,::Nothing) = nothing

@inline function unit_system(x1,x2,x3,x4)
    t1 = unit_system(unit_system(x1),unit_system(x2))
    t2 = unit_system(unit_system(x3),unit_system(x4))
    return unit_system(t1,t2)
end

unit_system(x::T,y::T) where T = x
unit_system(x::T,y::S) where {T,S}= throw(error("cannot use two unit systems at the same type."))

#fast escape for functions with no units
has_units(x) = Val(true)

#no units, default.
@inline with_output_unit(res,::Tuple{Nothing,Nothing},f::F) where F = res

#no units in the input, but output was provided
function with_output_unit(res,t::Tuple{Nothing,X},f::F) where {F,X}
    _,output = t
    return with_output_unit(has_units(f),res,(unit_system(output),output),f)
end

#no units in the output, but input has an unit system
function with_output_unit(res,t::Tuple{X,Nothing},f::F) where {F,X}
    unit_sys,_ = t
    return with_output_unit(has_units(f),res,(unit_sys,unit_type(unit_sys,f)),f)
end

with_output_unit(res::R,t::Tuple{X,Y},f::F) where {R,F,X,Y} = with_output_unit(Val(true),res,t,f)

with_output_unit(::Val{false},res,t,f) = res

#units provided by the input and by the output, check if we are in the same unit system.
function with_output_unit(::Val{true},res::R,t::Tuple{X,Y},f::F) where {R,F,X,Y}
    unit_sys0,output = t
    unit_sys = unit_system(unit_sys0,unit_system(output))
    output_from_f = unit_type(unit_sys,f)
    output_solved = solve_unit(unit_sys,output,output_from_f)
    return with_output_unit(res,(unit_sys,output_solved))
end

#fast escape for symbol properties
with_output_unit(::Val{true},res::Symbol,t::Tuple{X,Y},f::F) where {F,X,Y} = res

#differentiate p-T input from V-T input
unitful_is_pressure(p) = true