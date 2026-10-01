module DifferentialOperators;


export DiffOperator, Identity, Derivatex, Derivatey, AbstractGradient, Gradient, AdjointGradient, Divergence, Laplacian, DiffMatrix, AdjointDiffMatrix
export ∂x, ∂y, gradient, divergence, laplacian, ∇, Δ, diffmatrix


abstract type DiffOperator end
abstract type AbstractGradient <: DiffOperator end
abstract type AbstractIdentity <: DiffOperator end
abstract type AbstractDiffMatrix <: DiffOperator end
struct Identity <: AbstractIdentity end
# struct AdjointIdentity <: AbstractIdentity end
struct Derivatex <: DiffOperator end
struct Derivatey <: DiffOperator end
struct Gradient <: AbstractGradient end
struct AdjointGradient <: AbstractGradient end
struct Divergence <: DiffOperator end
struct Laplacian <: DiffOperator end
struct DiffMatrix <: AbstractDiffMatrix end
struct AdjointDiffMatrix <: AbstractDiffMatrix end

∂x = Derivatex()
∂y = Derivatey()
gradient = Gradient()
divergence = Divergence()
laplacian = Laplacian()
diffmatrix = DiffMatrix()
adjointdiffmatrix = AdjointDiffMatrix()


Base.adjoint(::Type{T}) where {T <: DiffOperator} = T
Base.adjoint(::Type{Gradient}) = AdjointGradient
Base.adjoint(::Type{AdjointGradient}) = Gradient
Base.adjoint(::Type{DiffMatrix}) = AdjointDiffMatrix
Base.adjoint(::Type{AdjointDiffMatrix}) = DiffMatrix

Base.adjoint(::Gradient) = AdjointGradient()
Base.adjoint(::AdjointGradient) = Gradient()
Base.adjoint(::DiffMatrix) = AdjointDiffMatrix()
Base.adjoint(::AdjointDiffMatrix) = DiffMatrix()

const ∇ = gradient
const Δ = laplacian

end
