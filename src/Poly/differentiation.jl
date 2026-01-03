
"""
```
   derivative(p::PolyField,z)
```

Compute the derivative of a PolyField with respect to the variable `z`.

# Examples
```
   julia> p = TensorPolynomial((1.,2.,3),(0.,2))
   (1.0 + 2.0*x + 3.0*x^2)(2.0*y)
   julia> pₓ = derivative(p,:x)
   (2.0 + 6.0*x)(2.0*y)
```
"""
function Polynomials.derivative(p::TensorPolynomial{F,X,Y},z::Symbol) where {F,X,Y}
    z == X && return TensorPolynomial(derivative(p.px),p.py)
    z == Y && return TensorPolynomial(p.px,derivative(p.py))
    throw(ArgumentError("Z must be an indeterminate present in the field, but Z=$z was provided for a field with indeterminates X=$X and Y=$Y"))
end

function Polynomials.derivative(p::PolySum{F,X,Y},z::Symbol) where {F,X,Y}
    z == X && return derivative(p.left,X) + derivative(p.right,X)
    z == Y && return derivative(p.left,X) + derivative(p.right,X)
    throw(ArgumentError("Z must be an indeterminate present in the field, but Z=$z was provided for a field with indeterminates X=$X and Y=$Y"))
end

Polynomials.derivative(p::PolyTensorField,z::Symbol) = PolyTensorField(derivative.(p.tensor,z))

"""
```
   gradient(p::PolyScalarField{F,X,Y}) where {F,X,Y}
```
Computes the gradient of a `PolyScalarField` and returns a `PolyVectorField`. 
"""
function gradient(p::PolyScalarField{F,X,Y}) where {F,X,Y}
    dx = derivative(p,X)
    dy = derivative(p,Y)
    PolyVectorField([dx,dy])
end

"""
```
   divergence(p::PolyVectorField{F,X,Y}) where {F,X,Y}
```
Computes the divergence of a `PolyVectorField` and returns a `PolyScalarField`, typically a  `PolySum`. 
"""
function divergence(v::PolyVectorField)
    d1x = derivative(v[1].px)
    d2y = derivative(v[2].py)
    part1 = TensorPolynomial(d1x,v[1].py)
    part2 = TensorPolynomial(v[2].px,d2y)
    PolySum(part1,part2)
end

LinearAlgebra.dot(::Type{gradient},v::PolyVectorField) = divergence(v)


"""
```
   laplacian(p::PolyScalarField{F,X,Y}) where {F,X,Y}
```
Computes the laplacian of a `PolyScalarField` and returns another `PolyScalarField`, typically a `PolySum`. 
"""
laplacian(v::PolyScalarField) = divergence(gradient(v))

# This method should be removed in the next Polynomials update
Polynomials.derivative(p::ImmutablePolynomial{F,X,1}) where {F,X} = zero(p)
