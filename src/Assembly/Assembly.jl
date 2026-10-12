module Assembly

using LinearAlgebra
using FixedSizeArrays
using Dictionaries
using SparseArrays
using Collects

using ..DifferentialOperators
using ..Meshes
using ..Tensors
using ..PolyFields
using ..Integration
using ..Forms
using ..Measures


DICT_OP = Dict([(*,⊗)])

include("reftensor.jl")
include("collapser.jl")
include("matrices.jl")

export RefTensor, Collapser
export assembly_matrix, assembly_rhs

end
