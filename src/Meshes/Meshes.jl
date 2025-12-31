module Meshes

using Dictionaries
using ExactPredicates
using LinearAlgebra
using Makie
using Markdown
using Printf
using Triangulate
using StaticArrays
using DocStringExtensions


export Edge,Triangle,HPMesh
export triangle,hpmesh
export inttype,floattype,degtype
export degree,degrees,psortperm,edges,triangles,psortednodes,psortededges
export tridof,edgedof
export circmesh,circmesh_graded_center,rectmesh,l_mesh,l_graded
export plothpmesh
export BoundaryHPMesh,dirichletboundary,neumannboundary
export set_dirichlet!, set_neumann!
export refine!
export mark!

BOUNDARY_DICT = Dict(:dirichlet=>1,:neumann=>2)

include("settuple.jl")
include("edge.jl")
include("triangle.jl")
include("mesh.jl")
include("boundary_mesh.jl")
include("refine.jl")
include("show.jl")
include("plots.jl")
include("some_meshes.jl")



const HASH_SEED = UInt === UInt64 ? 0x793bac59abf9a1da : 0xdea7f1da

end;