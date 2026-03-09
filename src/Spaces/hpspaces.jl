abstract type AbstractSpace end

abstract type ScalarSpace <: AbstractSpace end
abstract type VectorSpace <: AbstractSpace end
abstract type TensorSpace <: AbstractSpace end

struct StdScalarSpace <: ScalarSpace end
struct StdVectorSpace <: VectorSpace end

struct OperatorSpace{F<:DiffOperator, O<:Order, S <: AbstractSpace} <: AbstractSpace
    operator::F
    space::S
end

(::Gradient)(s::StdScalarSpace) = OperatorSpace(gradient, s)
(::Divergence)(s::Union{StdVectorSpace, OperatorSpace}) = OperatorSpace(divergence, s)
(::Laplacian)(s::StdScalarSpace) = OperatorSpace(laplacian, s)

basis(::StdScalarSpace, p::Tuple) = StandardBasis(p)

struct Operation{F, S1, S2}
    operator::F
    left::S1
    right::S2
end
Operation(a, b) = Operation{typeof(a), typeof(b), Nothing}(a, b, nothing)

const Sp = Union{AbstractSpace, Operation}

Base.:*(thing::Union{Number, AbstractArray, Function, PolyField}, s::Sp) = Operation(*, thing, s)
Base.:*(s::Sp, thing::Union{Number, AbstractArray, Function, PolyField}) = thing * s
Base.:*(s₁::Sp, s₂::Sp) = Operation(*, s₁, s₂)

LinearAlgebra.dot(a::AbstractArray, op::Sp) = Operation(dot, a, op)
LinearAlgebra.dot(op::Sp, a::AbstractArray) = Operation(dot, op, a)
LinearAlgebra.dot(a::Sp, b::Sp) = Operation(dot, a, b)

struct Order{B} end

order(::Any) = Order{0}()
order(::Union{AbstractSpace, Type{<:AbstractSpace}}) = Order{(0,)}()
order(::Divergence) = Order{1}()
order(::Gradient) = Order{1}()
order(::Laplacian) = Order{2}()

order(op::OperatorSpace) = order(order(op.operator), order(op.space))
function order(::Order{B}, ::Order{C}) where {B, C}
    if B isa Integer && C isa Integer
        Order{B + C}()
    elseif B isa Integer && C isa Tuple
        length(C) == 1 || ArgumentError("Combination not implemented")
        (D,) = C
        Order{(D + B,)}()
    else
        ArgumentError("Combination not implemented")
    end
end
order(op::Operation) = combine(order(op.left), order(op.right))
combine(::Order{B}, ::Order{C}) where {B, C} = Order{(B..., C...)}()
combine(::Order{B}, ::Order{0}) where {B} = Order{B}()
combine(::Order{0}, ::Order{B}) where {B} = Order{B}()
combine(::Order{0}, ::Order{0}) = Order{0}()


# order(::OperatorSpace{::typeof(gradient),::AbstractSpace})  = Order{(1,)}()
# # order(::OperatorSpace{::typeof(divergence),::AbstractSpace}) = Order{(1,)}()
# order(::OperatorSpace{::typeof(laplacian),::AbstractSpace})  = Order{(2,)}()
# function order(op::T) where T<:OperatorSpace{typeof(divergence)}
#     order(divergence,order(op.space))
# end
# order(::typeof(divergence),::Order{(B,)}) where B = Order{(B+1,)}()
# order(::OperatorSpace{typeof(divergence),::AbstractSpace},::Order{0}) = Order{(1,)}()
# function order(::OperatorSpace{::typeof(gradient),inner::OperatorSpace})
#     order(::typeof(gradient),order(inner))
# end
# order(::OperatorSpace{typeof(laplacian)},::Order{(B,)}) = Order{(B+1,)}()
# function order(::OperatorSpace{::typeof(laplacian)},inner::OperatorSpace})
#     order(::typeof(laplacian),order(inner))
# end
# order(::OperatorSpace{typeof(divergence)},::Order{(B,)}) = Order{(B+1,)}()


