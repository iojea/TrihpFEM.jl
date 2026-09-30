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

include("localtensor.jl")
include("matrices.jl")

export LocalTensor, _eval_operation, get_shape_functions, collapser

end
