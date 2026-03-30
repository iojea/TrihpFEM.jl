module Measures

using LinearAlgebra
using StaticArrays
using FixedSizeArrays
using Dictionaries

using ..Meshes
using ..PolyFields
using ..Integration

import ..Meshes: elements

include("auxiliarydata.jl")
include("meas.jl")

export Measure, elements

end; #module