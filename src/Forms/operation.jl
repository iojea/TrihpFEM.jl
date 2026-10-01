"""
    abstract type CoeffType end
An abstract type for defining a treat that allow `IntegrationTerm`s to be dispatched to different methods for integration. Exact integration for constant coefficient terms and quadrature rules for variable coefficient terms.
"""
abstract type CoeffType end

struct ConstantCoeff <: CoeffType end
struct VariableCoeff <: CoeffType end

coefftype(_) = ConstantCoeff
coefftype(::Function) = VariableCoeff
Base.promote_type(::Type{ConstantCoeff}, ::Type{VariableCoeff}) = VariableCoeff


struct Operation{C <: CoeffType, U <: Tuple}
    parts::U
    Operation{C}(t...) where {C <: CoeffType} = new{C, typeof(t)}(t)
end

Base.getindex(op::Operation, i) = getindex(op.parts, i)
Base.length(op::Operation) = length(op.parts)
Base.lastindex(op::Operation) = length(op)
const NumArray = Union{Number, AbstractArray}

function unroll(sf::ShapeFunction)
    return isadjoint(sf) ? (sf, *, chain_operator(sf)) : (chain_operator(sf), *, sf)
end


Base.:*(a::NumArray, b::ShapeFunction) = Operation{ConstantCoeff}(a, *, unroll(b)...)
Base.:*(a::ShapeFunction, b::NumArray) = Operation{ConstantCoeff}(unroll(b)..., *, a)
Base.:*(a::ShapeFunction, b::ShapeFunction) = Operation{ConstantCoeff}(unroll(a)..., *, unroll(b)...)
Base.:*(a::Function, b::ShapeFunction) = Operation{VariableCoeff}(a, *, unroll(b)...)
Base.:*(a::ShapeFunction, b::Function) = Operation{VariableCoeff}(unroll(a)..., *, b)
Base.:*(a::Operation{C}, b::ShapeFunction) where {C <: CoeffType} = Operation{C}(a.parts..., *, unroll(b)...)
Base.:*(a::ShapeFunction, b::Operation{C}) where {C <: CoeffType} = Operation{C}(unroll(a)..., *, b)
Base.:*(a::NumArray, b::Operation{C}) where {C <: CoeffType} = Operation{C}(a, *, b.parts...)
Base.:*(a::Operation{C}, b::NumArray) where {C <: CoeffType} = Operation{C}(a.parts..., *, b)
function Base.:*(a::Operation{C₁}, b::Operation{C₂}) where {C₁ <: CoeffType, C₂ <: CoeffType}
    C = promote_type(C₁, C₂)
    return Operation{C}(a.parts..., *, b.parts...)
end

_adj(a::Any) = a'
_adj(f::Function) = adjoint ∘ f
_adj(::typeof(*)) = *
_adj(::typeof(LinearAlgebra.dot)) = LinearAlgebra.dot

Base.adjoint(op::Operation{C}) where {C <: CoeffType} = Operation{C}(reverse(_adj.(op.parts))...)

LinearAlgebra.dot(a::NumArray, b::ShapeFunction) = a' * b
LinearAlgebra.dot(a::ShapeFunction, b::NumArray) = a' * b
LinearAlgebra.dot(a::ShapeFunction, b::ShapeFunction) = a' * b
LinearAlgebra.dot(a::Function, b::ShapeFunction) = (transpose ∘ a) * b
LinearAlgebra.dot(a::ShapeFunction, b::Function) = a' * b
LinearAlgebra.dot(a::Operation, b::ShapeFunction) = a' * b
LinearAlgebra.dot(a::ShapeFunction, b::Operation) = a' * b
LinearAlgebra.dot(a::Operation, b::Operation) = a' * b

Base.:-(u::ShapeFunction) = (-1) * u


"""
   Integrand{C<:CoeffType,T<Tuple,N}
An expression for integration, depending on at least one ShapeFunction. `N` is the number of `ShapeFunction`s involved  `C` the type of coefficient of the factor. `T` is a type of the form `Tuple{O1}` o `Tuple{O1,O2}`, where `O1,O2<:DiffOperator`, indicates the differential operators applied to the `ShapeFunction`s involved."""
struct Integrand{C <: CoeffType, T <: Tuple, N}
    operation::Operation{C, T}
    function Integrand(operation::Operation{C, T}) where {C <: CoeffType, T <: Tuple}
        return new{C, T, length(T.parameters)}(operation)
    end
end
coefftype(::Integrand{C}) where {C} = C()


Integrand(u::ShapeFunction) = 1 * u

const ∫ = Integrand
