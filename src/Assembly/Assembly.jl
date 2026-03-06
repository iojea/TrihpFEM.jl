module Assembly

    using LinearAlgebra
    using FixedSizeArrays
    using EllipsisNotation
    using Dictionaries
    using SparseArrays
    using Tensors

    using ..Meshes
    using ..Poly
    using ..Spaces
    using ..Integration
    using ..Forms
    using ..Measures

    include("matrices.jl")

    export integrate
end; #module
