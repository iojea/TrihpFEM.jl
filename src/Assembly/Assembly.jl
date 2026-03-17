module Assembly

using LinearAlgebra
using FixedSizeArrays
using EllipsisNotation
using Dictionaries
using SparseArrays
using Tensors
using Collects

using ..Meshes
using ..PolyFields
using ..Integration
using ..Forms
using ..Measures

include("matrices.jl")

export ref_tensors, assembly_matrix

end
