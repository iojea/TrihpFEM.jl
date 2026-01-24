module Measures

using LinearAlgebra
using StaticArrays
using Dictionaries
using ..Meshes
using ..Integration

include("auxiliarydata.jl")
include("meas.jl")

export Measure

end; #module