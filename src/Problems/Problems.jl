module Problems

using CommonSolve
using StaticArrays
using Dictionaries
using FixedSizeArrays
using LinearAlgebra
using ..Meshes
using ..Poly
using ..Spaces
using ..Integration
using ..Forms


include("auxiliarydata.jl")
include("measures.jl")
include("feproblem.jl")


export Measure
export FEProblem

end; #module
