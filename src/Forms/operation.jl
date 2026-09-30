
"""
    abstract type CoeffType end
An abstract type for defining a treat that allow `IntegrationTerm`s to be dispatched to different methods for integration. Exact integration for constant coefficient terms and quadrature rules for variable coefficient terms.
"""
abstract type CoeffType end

struct ConstantCoeff <: CoeffType end
struct VariableCoeff <: CoeffType end

coefftype(_) = ConstantCoeff
coefftype(::Function) = VariableCoeff
Base.promote_type(::Type{ConstantCoeff},::Type{VariableCoeff}) = VariableCoeff



struct Operation{C<:CoeffType,T<:Tuple,OP<:Function,U<:Tuple}
    op::OP
    parts::U
    function Operation{C,T}(op,a,b) where {C<:CoeffType,T<:Tuple}
        any(hasoperator.((a,b))) || error("A `ShapeFunction` must be present in an `Operation`")
        new{C,T,typeof(op),typeof((a,b))}(op,(a,b))
    end  
end

hasoperator(::Any) = false
hasoperator(::ShapeFunction) = true
hasoperator(::Operation) = true


function Base.:*(a::A,b::B) where {A<:Union{Number,AbstractArray},O,B<:ShapeFunction{O}}
    Operation{ConstantCoeff,Tuple{O}}(*,a,b)
end
function Base.:*(a::A,b::B) where {O,A<:ShapeFunction{O},B<:Union{Number,AbstractArray}}
    Operation{ConstantCoeff,Tuple{O}}(*,a,b)
end
function Base.:*(a::ShapeFunction{O1},b::ShapeFunction{O2}) where {O1,O2}
    Operation{ConstantCoeff,Tuple{O1,O2}}(*,a,b)
end
function Base.:*(a::A,b::B) where {A<:Function,O,B<:ShapeFunction{O}}
    Operation{VariableCoeff,O}(*,a,b)
end
function Base.:*(a::A,b::B) where {O,A<:ShapeFunction{O},B<:Function}
    Operation{VariableCoeff,O}(*,a,b)
end
function Base.:*(a::Operation{C,Tuple{O1}},b::ShapeFunction{O2}) where {C,O1,O2}
    Operation{C,Tuple{O1,O2}}(*,a,b)
end
function Base.:*(a::ShapeFunction{O1},b::Operation{C,Tuple{O2}}) where {C,O1,O2}
    Operation{C,Tuple{O1,O2}}(*,a,b)
end
function Base.:*(a::Union{Number,AbstractArray},b::Operation{C,T}) where {C,T}
    Operation{C,T}(*,a,b)
end
function Base.:*(b::Operation{C,T},a::Union{Number,AbstractArray}) where {C,T}
    Operation{C,T}(*,b,a)
end
function Base.:*(a::Operation{C₁,Tuple{T₁}},b::Operation{C₂,Tuple{T₂}}) where {C₁,C₂,T₁,T₂}
    C = promote_type(C₁,C₂)
    Operation{C,Tuple{T₁,T₂}}(*,a,b)
end

_adj(a::Any) = a'
_adj(f::Function) = adjoint∘f

Base.adjoint(op::Operation{C,T,typeof(*)}) where {C,T} = _adj(op.parts[2])*_adj(op.parts[1])

LinearAlgebra.dot(a::A,b::B) where {A<:Union{Number,AbstractArray},B<:ShapeFunction} = a'*b
LinearAlgebra.dot(a::A,b::B) where {A<:ShapeFunction,B<:Union{Number,AbstractArray}} = a'*b
LinearAlgebra.dot(a::ShapeFunction,b::ShapeFunction) = a'*b
LinearAlgebra.dot(a::A,b::B) where {A<:Function,B<:ShapeFunction} = (transpose∘a)*b
LinearAlgebra.dot(a::A,b::B) where {A<:ShapeFunction,B<:Function} = a'*b
LinearAlgebra.dot(a::Operation,b::ShapeFunction) = a'*b
LinearAlgebra.dot(a::ShapeFunction,b::Operation) = a'*b
LinearAlgebra.dot(a::Operation,b::Operation) = a'*b

Base.:-(u::ShapeFunction) = (-1)*u

"""
   Integrand{C<:CoeffType,T<Tuple,N}
An expression for integration, depending on at least one ShapeFunction. `N` is the number of `ShapeFunction`s involved  `C` the type of coefficient of the factor. `T` is a type of the form `Tuple{O1}` o `Tuple{O1,O2}`, where `O1,O2<:DiffOperator`, indicates the differential operators applied to the `ShapeFunction`s involved."""
struct Integrand{C<:CoeffType,T<:Tuple,N}
    operation::Operation{C,T}
    function Integrand(operation::Operation{C,T}) where {C<:CoeffType,T<:Tuple}
        new{C,T,length(T.parameters)}(operation)
    end
end
coefftype(::Integrand{C}) where C = C()


Integrand(u::ShapeFunction) = 1*u

const ∫ = Integrand


