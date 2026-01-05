abstract type AbstractSpace end

abstract type ScalarSpace <: AbstractSpace end
abstract type VectorSpace <: AbstractSpace end
abstract type TensorSpace <: AbstractSpace end

struct StdScalarSpace <: ScalarSpace end
struct StdVectorSpace <: VectorSpace end
struct StdTensorSpace <: TensorSpace end

struct OperatorSpace{F<:Function,S<:AbstractSpace} <: AbstractSpace
    operator::F
    space::S
end


Poly.gradient(s::StdScalarSpace) = OperatorSpace(gradient,s)
Poly.divergence(s::StdVectorSpace) = OperatorSpace(divergence,s)
Poly.laplacian(s::StdScalarSpace) = OperatorSpace(laplacian,s)

struct Order{B} end

order(_) = Order{0}()
order(::AbstractSpace) = Order{(0,)}()
order(::OperatorSpace{typeof(gradient),S}) where S = Order{(1,)}()
order(::OperatorSpace{typeof(divergence),S}) where S = Order{(1,)}()
order(::OperatorSpace{typeof(laplacian),S}) where S = Order{(2,)}()
combine(::Order{B},::Order{C}) where {B,C} = Order{(B...,C...)}()
combine(::Order{B},::Order{0}) where B = Order{B}()
combine(::Order{0},::Order{B}) where B = Order{B}()


basis(::StdScalarSpace,p::Tuple) = StandardBasis(p)

# A trait for identifying Constant Coefficients, which allow precomputation of local tensors.
abstract type CoeffType end

struct Constant <: CoeffType end
struct Variable <: CoeffType end

coefftype(::AbstractSpace) = Constant()
coefftype(::AbstractArray) = Constant()
coefftype(::Number) = Constant()
coefftype(::Function) = Variable()

Base.promote(::Constant,::Variable) = Variable()
Base.promote(::Variable,::Constant) = Variable()
Base.promote(::Variable,::Variable) = Variable()
Base.promote(::Constant,::Constant) = Constant()
Base.promote(::Nothing,::B) where B<: CoeffType = B()
Base.promote(::B,::Nothing) where B<: CoeffType = B()

Base.promote_type(::Type{Constant},::Type{Variable}) = Variable
Base.promote_type(::Type{Variable},::Type{Variable}) = Variable
Base.promote_type(::Type{Constant},::Type{Constant}) = Constant
Base.promote_type(::Type{Variable},::Type{Constant}) = Variable
Base.promote_type(::Type{Nothing},::Type{B}) where B<:CoeffType = B
Base.promote_type(::Type{B},::Type{Nothing}) where B<:CoeffType = B


struct Operation{F<:Function,S1,S2}
    operator::F
    left::S1
    right::S2
end
Operation(a,c,b)  = Operation{typeof(a),typeof(b),typeof(c)}(a,b,c)
Operation(a,b) = Operation{typeof(a),typeof(b),Nothing}(a,b,nothing)

const Sp = Union{AbstractSpace,Operation}

Base.:*(thing::Union{Number,AbstractArray,Function,PolyField},s::Sp) = Operation(*,thing,s)
Base.:*(s::Sp,thing::Union{Number,AbstractArray,Function,PolyField}) = thing*s
Base.:*(s₁::Sp,s₂::Sp) = Operation(*,s₁,s₂)

LinearAlgebra.dot(a::AbstractArray,op::Sp) = Operation(dot,a,op)
LinearAlgebra.dot(op::Sp,a::AbstractArray) = Operation(dot,op,a)
LinearAlgebra.dot(a::Sp,b::Sp) = Operation(dot,a,b)

coefftype(op::Operation) = promote(coefftype(op.left),coefftype(op.right))

order(op::Operation) = combine(order(op.left),order(op.right))


#Integrand(op) = Integrand()

# Base.:*(f::Integran,m::Measur) = integrate(CoeffType(f.op),f.op,


# get_space(space::AbstractHPSpace) = space
# function get_space(opfield::OperationField)
#     (;args) = opfield
#     if typeof(args[1])<:AbstractHPSpace
#         return args[1]
#     elseif typeof(args[2])<:AbstractHPSpace
#         return args[2]
#     else
#         error("Something is wrong. An OperatorField was found with no space involved")
#     end
# end
# get_spaces(itd::Integrand) = get_space.(itd.func.args)

# function integrate(::Constant,f::Integrand,m::Measure)
#     M = build_tensors(f,m)
# end


# #function build_tensors()