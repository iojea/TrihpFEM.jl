module PolyFields

    using StaticArrays
    using Polynomials
    using LinearAlgebra
    using FixedSizeArrays
    using ..Meshes
    using ..Tensors

    import ..Tensors: otimes
    import ..DifferentialOperators: DiffOperator, Identity,Derivatex,Derivatey,AbstractGradient,Gradient,AdjointGradient, Divergence,Laplacian,DiffMatrix,AdjointDiffMatrix
    import ..DifferentialOperators: ∂x,∂y,gradient,∇,divergence,laplacian,Δ,diffmatrix,adjointdiffmatrix

    const VARIABLE_NAMES = (:x,:y,:z)
    const VARIABLE_IDX = Dict(:x=>1,:y=>2,:z=>3)
    
    include("fields.jl")
    include("polytensors.jl")
    include("affine.jl")
    include("differentiation.jl")
    include("legendre.jl")
    include("show.jl")

    export AbstractField
    export ProductPoly
    export PolyScalarField, PolyTensorField, PolyVectorField, PolyMatrixField
    export PolySum
    export indeterminate, indeterminates, degs
    export LegendreIterator, StandardBasis
    export dot
    export AffineToRef
    # export affine!
    export jac
    export area
    export EvalType, Eval, Compose, Pass
    export evaluate
    export otimes

end; #module
