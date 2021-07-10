
"""
nested_eltype(x)

returns the element type of a variable
on numbers, returns `typeof(x)`
on arrays, works recursively (first element) until reaching an array of numbers, and then returns `eltype(x)`
"""
function nested_eltype end

nested_eltype(x::Number) = typeof(x)
nested_eltype(x::AbstractArray{<:Number}) = eltype(x)
nested_eltype(x::AbstractArray) = nested_eltype(first(x))


"""
nested_copyto!(dest,src)

assigns the contents of `src` on `dest`, works on deeply nested arrays

"""
function nested_copyto! end

function nested_copyto!(dest::Number,src::Number)
return dest 
end

function nested_copyto!(dest::T,src::T) where T<:AbstractArray{<:Number}
copyto!(dest,src)
return dest
end

function nested_copyto!(dest::T,src::T) where T<:AbstractArray
nested_copyto!.(dest,src)
return dest
end


    