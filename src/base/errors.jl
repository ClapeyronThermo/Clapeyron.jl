not_eos_error(model) = throw(ArgumentError("$model does not have eos defined "))

"""
    single_component_check(method,model)

Checks if a model is a single component model, throws an error otherwise.
"""
function single_component_check(method,model)
    l = length(model)
    l == 1 && return nothing
    single_component_error(method,model)
end

function single_component_error(method,model)
    l = length(model)
    msg = string(method," only supports single component models, ",model," has ",l," components.")
    throw(DimensionMismatch(msg))
end

"""
    binary_component_check(method,model)

Checks if a model is a single component model, throws an error otherwise.
"""
function binary_component_check(method,model)
    l = length(model)
    l == 2 && return nothing
    binary_component_error(method,model)
end

function binary_component_error(method,model)
    l = length(model)
    msg = string(method," only supports binary component models, ",model," has ",l," components.")
    throw(DimensionMismatch(msg))
end

function incorrect_squarematrix_error(model,n)
    l = length(model)
    msg = string(model," has $l components, while input matrix is of size $(n)×$(n)")
    throw(DimensionMismatch(msg))
end
