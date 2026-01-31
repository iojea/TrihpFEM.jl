module Spaces

    using StaticArrays
    using Dictionaries
    using LinearAlgebra
    using Polynomials
    using ..Meshes
    using ..Poly
    include("hpspaces.jl")


    export AbstractSpace
    export StdScalarSpace, StdVectorSpace
    export ScalarSpace, VectorSpace
    export basis
    export OperatorSpace
    export gradient, divergence, laplacian
    export ∇, Δ
    export dot
    export Order
    export order


end; #module
