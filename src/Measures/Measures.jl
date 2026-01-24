module Measures

using LinearAlgebra
using StaticArrays
using FixedSizeArrays
using Dictionaries

using ..Meshes
using ..Poly
using ..Integration

include("auxiliarydata.jl")
include("meas.jl")

export Measure

end; #module