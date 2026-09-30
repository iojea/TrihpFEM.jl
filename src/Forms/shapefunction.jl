"""
    ShapeFunctionKind
A trait for differentiating Trial and Test ShapeFunctions. 
"""
abstract type ShapeFunctionKind end
struct Trial <: ShapeFunctionKind end
struct Test <: ShapeFunctionKind end


"""
    ShapeFunction(dim)
Defines a standard shape function, of dimension `dim`.
```
julia> u = ShapeFunction(1)
```
A gradient shape function can be build with `ShapeFunction{Gradient,1}()`. This means: the gradient of a shape function with `dim=1`. However the preferred way to build `ShapeFunction`s of this kind is to a apply the differential operator to the standard `ShapeFunction`:
```
julia> grad_u = ∇(u) 
```
"""
struct ShapeFunction{T <: ShapeFunctionKind, O <: DiffOperator, N, D}
    function ShapeFunction{T, O, N, D}() where {T <: ShapeFunctionKind, O <: DiffOperator, N, D}
        return new{T, O, N, D}()
    end
    function ShapeFunction{T}(; N = 1, D = 2) where {T <: ShapeFunctionKind}
        (D isa Integer && 1 <= D <= 2) || throw(ArgumentError("Only `ShapeFunction`s of dimension `1` or `2` can be built."))
        N isa Integer || throw(ArgumentError("Only integer dimensions can be used for defining a `ShapeFunction`. N=$N was passed."))
        N > 1 && throw(ArgumentError("Multidimensional `ShapeFunction`s are not implemented yet. Stay tuned."))
        return new{T, Identity, N, D}()
    end
end
ShapeFunction(kind::ShapeFunctionKind; N = 1, D = 2) = ShapeFunction{typeof(kind)}(; N = N, D = D)


function (::Gradient)(::ShapeFunction{T, Identity, 1, D}) where {T <: ShapeFunctionKind, D}
    return ShapeFunction{T, Gradient, 1, D}()
end
function (::Laplacian)(::ShapeFunction{T, Identity, 1, D}) where {T <: ShapeFunctionKind, D}
    return ShapeFunction{T, Laplacian, 1, D}()
end
function (::Derivatex)(::ShapeFunction{T, Identity, 1, D}) where {T <: ShapeFunctionKind, D}
    return ShapeFunction{T, Derivatex, 1, D}()
end
function (::Derivatey)(::ShapeFunction{T, Identity, 1, D}) where {T <: ShapeFunctionKind, D}
    return ShapeFunction{T, Derivatey, 1, D}()
end
(::Divergence)(::ShapeFunction{T, Gradient, 2, D}) where {T, D} = ShapeFunction{T, Laplacian, 1, D}()

operator(::ShapeFunction{T, O}) where {T, O} = O()
operatortype(::ShapeFunction{T, O}) where {T, O} = O

chain_operator(::ShapeFunction{T, Gradient, 1, D}) where {T, D} = DiffMatrix()
chain_operator(::ShapeFunction{T, AdjointGradient, 1, D}) where {T, D} = AdjointDiffMatrix()
chain_operator(::ShapeFunction{T, Identity, 1, D}) where {T, D} = Identity()

isadjoint(sf::ShapeFunction) = isadjoint(operator(sf))
isadjoint(::DiffOperator) = false
isadjoint(::AdjointGradient) = true
isadjoint(::AdjointDiffMatrix) = true

dim(::ShapeFunction{T, O, N, D}) where {T, O, N, D} = D
basis(::ShapeFunction{T, O, 1, 2}, degs) where {T, O} = StandardBasis(degs)
basis(::ShapeFunction{T, O, 1, 1}, deg) where {T, O} = LegendreIterator(deg)

Base.adjoint(::ShapeFunction{T, O, N, D}) where {T, O, N, D} = ShapeFunction{T, adjoint(O), N, D}()
