abstract type DiffOperator end

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


"""
```
   derivative(p::PolyField,z)
```

Compute the derivative of a PolyField with respect to the variable `z`.

# Examples
```
   julia> p = BiPoly((1.,2.,3),(0.,2))
   (1.0 + 2.0*x + 3.0*x^2)(2.0*y)
   julia> pₓ = derivative(p,:x)
   (2.0 + 6.0*x)(2.0*y)
```
"""
function Polynomials.derivative(p::BiPoly{F, X, Y}, z::Symbol) where {F, X, Y}
    z == X && return BiPoly(derivative(p.px), p.py, X, Y)
    z == Y && return BiPoly(p.px, derivative(p.py), X, Y)
    throw(ArgumentError("$z is not an indeterminate of the polynomial"))
end

function Polynomials.derivative(p::PolySum{F, X, Y}, z::Symbol) where {F, X, Y}
    z == X && return derivative(p.left, X) + derivative(p.right, X)
    z == Y && return derivative(p.left, Y) + derivative(p.right, Y)
    throw(ArgumentError("$z is not an indeterminate of the polynomial"))
end

Polynomials.derivative(p::PolyTensorField, z::Symbol) = PolyTensorField(derivative.(p.tensor, z))

(::Derivatex)(p::PolyField{F,X,Y}) where {F,X,Y} = derivative(p,X)
(::Derivatey)(p::PolyField{F,X,Y}) where {F,X,Y} = derivative(p,Y)
              
"""
```
   gradient(p::PolyScalarField{F,X,Y}) where {F,X,Y}
```
Computes the gradient of a `PolyScalarField` and returns a `PolyVectorField`. 
"""
(::Gradient)(p::PolyScalarField) = PolyVectorField([∂x(p),∂y(p)])

"""
```
   divergence(p::PolyVectorField{F,X,Y}) where {F,X,Y}
```
Computes the divergence of a `PolyVectorField` and returns a `PolyScalarField`, typically a  `PolySum`. 
"""
(::Divergence)(v::PolyVectorField) = ∂x(v[1])+∂y(v[2])


"""
```
   laplacian(p::PolyScalarField{F,X,Y}) where {F,X,Y}
```
Computes the laplacian of a `PolyScalarField` and returns another `PolyScalarField`, typically a `PolySum`. 
"""
laplacian(v::PolyScalarField) = divergence(gradient(v))

# This method should be removed in the next Polynomials update
Polynomials.derivative(p::ImmutablePolynomial{F, X, 1}) where {F, X} = zero(p)
