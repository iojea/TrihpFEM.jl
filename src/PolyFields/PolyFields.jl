module PolyFields

    using StaticArrays
    using Polynomials
    using LinearAlgebra
    using FixedSizeArrays
    using Tensors
    using ..Meshes

    import ..DifferentialOperators: DiffOperator, Identity,Derivatex,Derivatey,Gradient,Divergence,Laplacian
    import ..DifferentialOperators: ∂x,∂y,gradient,∇,divergence,laplacian,Δ

    include("fields.jl")
    include("differentiation.jl")
    include("legendre.jl")
    include("show.jl")

    export PolyField
    export BiPoly
    export PolyScalarField, PolyVectorField, PolyTensorField
    export GeneralField
    export PolySum
    export indeterminate, indeterminates, degs
    export LegendreIterator, StandardBasis
    export dot
    export AffineToRef
    export affine!
    export jac
    export area
    export EvalType, Eval, Compose, Pass
    export evaluate

end; #module
