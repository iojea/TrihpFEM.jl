module Assembly

using LinearAlgebra
using TensorOperations
using FixedSizeArrays
using EllipsisNotation
using Dictionaries
using SparseArrays


using ..Meshes
using ..Poly
using ..Spaces
using ..Integration
using ..Forms
using ..Measures

include("matrices.jl")

export integrate
end; #module