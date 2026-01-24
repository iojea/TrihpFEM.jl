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

"""
```
   gradient(p::PolyScalarField{F,X,Y}) where {F,X,Y}
```
Computes the gradient of a `PolyScalarField` and returns a `PolyVectorField`. 
"""
function gradient(p::PolyScalarField{F, X, Y}) where {F, X, Y}
    dx = derivative(p, X)
    dy = derivative(p, Y)
    return PolyVectorField([dx, dy])
end

"""
```
   divergence(p::PolyVectorField{F,X,Y}) where {F,X,Y}
```
Computes the divergence of a `PolyVectorField` and returns a `PolyScalarField`, typically a  `PolySum`. 
"""
function divergence(v::PolyVectorField)
    X, Y = indeterminates(v)
    d1x = derivative(v[1], :x)
    d2y = derivative(v[2], :y)
    return d1x + d2y
end


"""
```
   laplacian(p::PolyScalarField{F,X,Y}) where {F,X,Y}
```
Computes the laplacian of a `PolyScalarField` and returns another `PolyScalarField`, typically a `PolySum`. 
"""
laplacian(v::PolyScalarField) = divergence(gradient(v))

# This method should be removed in the next Polynomials update
Polynomials.derivative(p::ImmutablePolynomial{F, X, 1}) where {F, X} = zero(p)
