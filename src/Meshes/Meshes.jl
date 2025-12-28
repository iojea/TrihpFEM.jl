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


export Edge,Triangle

# export degrees,psortperm,edges,pnodes,pedges,triangles,edges,degree
# export meshhp,circmesh,circmesh_graded_center,rectmesh
# export show, plothpmesh, animate_refinement
# export domainmesh
# export degrees_of_freedom!
# export inttype, floattype, degtype
# export tridofs


include("settuple.jl")
include("edge.jl")
include("triangle.jl")
# include("mesh.jl")
# include("boundary_mesh.jl")
# include("refine.jl")
# include("show.jl")
# include("plots.jl")
# include("examples.jl")

const HASH_SEED = UInt === UInt64 ? 0x793bac59abf9a1da : 0xdea7f1da

end;