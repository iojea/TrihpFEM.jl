(::Identity)(p::PolyScalarField) = p

"""
```
   derivative(p::AbstractField,z)
```

Compute the derivative of a AbstractField with respect to the variable `z`.

# Examples
```
   julia> p = ProductPoly((1.,2.,3),(0.,2))
   (1.0 + 2.0*x + 3.0*x^2)(2.0*y)
   julia> pₓ = derivative(p,:x)
   (2.0 + 6.0*x)(2.0*y)
```
"""
Polynomials.derivative(p::ProductPoly,::Val{:x}) = ProductPoly(derivative(p.polys[1]),p.polys[2:end]...)
Polynomials.derivative(p::ProductPoly{2},::Val{:y}) = ProductPoly(p.polys[1],derivative(p.polys[2]))
Polynomials.derivative(p::ProductPoly{3},::Val{:y}) = ProductPoly(p.polys[1],derivative(p.polys[2]),p.polys[3])
Polynomials.derivative(p::ProductPoly{3},::Val{:z}) = ProductPoly(p.polys[1:2]...,derivative(p.polys[3]))
Polynomials.derivative(p::ProductPoly,s::Symbol) = derivative(p,Val(s))

Polynomials.derivative(p::PolySum,s::Symbol) = derivative(p.left,s)+derivative(p.right,s)

(::Derivatex)(p::PolyScalarField) = derivative(p, :x)
(::Derivatey)(p::PolyScalarField) = derivative(p, :y)

"""
```
   gradient(p::PolyScalarField{F,X,Y}) where {F,X,Y}
```
Computes the gradient of a `PolyScalarField` and returns a `PolyVectorField`. 
"""
(::AbstractGradient)(p::PolyScalarField{2,F}) where {F} = Tensor{1,2}(∂x(p), ∂y(p))

"""
```
   divergence(p::PolyVectorField{F,X,Y}) where {F,X,Y}
```
Computes the divergence of a `PolyVectorField` and returns a `PolyScalarField`, typically a  `PolySum`. 
"""
(::Divergence)(v::PolyVectorField) = ∂x(v[1]) + ∂y(v[2])


"""
```
   laplacian(p::PolyScalarField{F,X,Y}) where {F,X,Y}
```
Computes the laplacian of a `PolyScalarField` and returns another `PolyScalarField`, typically a `PolySum`. 
"""
(::Laplacian)(v::PolyScalarField) = divergence(gradient(v))


function (::DiffMatrix)(v::Tensor{1,2,T}) where {T<:PolyScalarField}
   Tensor{2,2}(∂x(v[1]),∂x(v[2]),∂y(v[1]),∂y(v[2]))
end

(::DiffMatrix)(aff::AffineToRef) = inv(aff.A)
(::AdjointDiffMatrix)(aff::AffineToRef) = inv(aff.A)'
(::Identity)(::AffineToRef{D,F}) where {D,F} = one(F)

