module DifferentialOperators;

export DiffOperator, Identity, Derivatex, Derivatey, Gradient, Divergence, Laplacian
export ∂x,∂y,gradient,divergence,laplacian,∇,Δ


abstract type DiffOperator end

struct Identity <: DiffOperator end
struct Derivatex <: DiffOperator end
struct Derivatey <: DiffOperator end
struct Gradient <: DiffOperator end
struct Divergence <: DiffOperator end
struct Laplacian <: DiffOperator end


∂x = Derivatex()
∂y = Derivatey()
gradient = Gradient()
divergence = Divergence()
laplacian = Laplacian()

const ∇ = gradient
const Δ = laplacian

end