module Assembly

using LinearAlgebra
using FixedSizeArrays
using EllipsisNotation
using Dictionaries
using SparseArrays
using Tensors
using Collects

using ..DifferentialOperators
using ..Meshes
using ..PolyFields
using ..Integration
using ..Forms
using ..Measures

include("matrices.jl")

export assembly_matrix, assembly_rhs

end
