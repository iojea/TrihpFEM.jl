module Meshes

    using StyledStrings
    using Dictionaries
    using ExactPredicates
    using LinearAlgebra
    using Makie
    using Markdown
    using Printf
    using Triangulate
    using StaticArrays
    using DocStringExtensions


    export Point,Edge, EdgeAttributes, Triangle, TriangleAttributes, HPTriangulation, HPMesh,BoundaryHPMesh, DOF
    export triangle, hpmesh, data, elements
    export tag, degree, dof, ndof, longestedge, tagged_dof
    export ismarked, isgreen, isblue, isred, istagged, isinterior, isboundary
    export inttype, floattype, degtype
    export degrees_of_freedom!, isempty, empty!
    export settag!, setdegree!, setboundary!, setdirichlet!, setneumann!
    export degrees, psortperm, edges, triangles, psortednodes, psortededges
    export plothpmesh
    export BoundaryHPMesh, dirichletboundary, neumannboundary, domainmesh
    export mark!, refine!, p_conformity!, check_p_conformity, setdegrees!
    export circmesh, circmesh_graded_center, rectmesh, squaremesh
    export plothpmesh, plothpmesh2, degplot
    export details


    const BOUNDARY_DICT = Dict(:dirichlet => 1, :neumann => 2)
    const COLOR_DICT = Dict(0=>:lightgrey,1=>:forestgreen,2=>:cornflowerblue,3=>:brown3)
    const FACE_DICT = Dict(
                0=>StyledStrings.Face(foreground=StyledStrings.SimpleColor(212,212,212)),
                1=>StyledStrings.Face(foreground=StyledStrings.SimpleColor(33,140,33)),
                2=>StyledStrings.Face(foreground=StyledStrings.SimpleColor(99,148,237)),
                3=>StyledStrings.Face(foreground=StyledStrings.SimpleColor(204,51,51)),
                                                  )
    
    const TAG_DICT = Dict(0=>"Ω°",1=>"∂𝔇",2=>"∂𝔑")
    const _EDGE_TAG_COLORS = Dict(0 => :white, 1 => :cornflowerblue, 2 => :seagreen, 3 => :orange)
    const _MARKED_EDGE_COLOR = :black

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
