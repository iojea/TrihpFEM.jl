    using ..TrihpFEM.Meshes
    using Test
    using Dictionaries
    using StaticArrays

    #Edge creation and comparison.
    e₁ = Edge(UInt8(3),4)
    e₂ = Edge((4,3))
    e₃ = Edge(SVector(4,3))
    e₄ = Edge{Int32}(3,4)
    e₅ = Edge(e for e in e₄)
    @test repr(e₁) == "(3, 4)"
    @test repr(e₂) == "(4, 3)"
    @test isequal(e₁,e₂) && isequal(e₁,e₃) && isequal(e₄,e₅)
    @test typeof(e₄) == Edge{Int32}
    
    ea₁ = EdgeAttributes{UInt8}(1,0,false)
    @test repr(ea₁) == "(0x01, :Ω°, :noref)"
    de = Dictionary([e₁],[ea₁])
    @test e₁ in keys(de)
    @test de[e₂] == ea₁
    @test degree(de[e₁]) == 1
    @test !Meshes.ismarked(de[e₁])
    @test tag(de[e₁]) == 0
    mark!(de[e₁])
    @test Meshes.ismarked(de[e₂])
    setdegree!(de[e₁],3)
    @test typeof(tag(de[e₁])) == UInt8
    @test degree(de[e₁]) == 3
    @test isinterior(de[e₂])
    @test data(e₁) == (3,4)

    #Triangle creation and comparison
    t₁ = Triangle(Int32(1),2,3)
    t₂ = Triangle((2,3,1))
    t₃ = Triangle(SVector(3,1,2))
    t₄ = Triangle{UInt32}(2,1,3)
    t₅ = Triangle(t for t in t₄)
    @test isequal(t₁,t₂) && isequal(t₁,t₃) && isequal(t₁,t₄) && isequal(t₄,t₅)
    @test eltype(data(t₄)) == UInt32

    ta₁ = TriangleAttributes{UInt8,Float64}()
    ta₂ = TriangleAttributes()
    @test all(getproperty(ta₁,prop)[]==getproperty(ta₂,prop)[] for prop in propertynames(ta₁))

    dt = Dictionary([t₁],[ta₁])
    @test t₂ in keys(dt)
    @test !Meshes.ismarked(dt[t₂])
    mark!(dt[t₁],1)
    @test Meshes.ismarked(dt[t₂])
    @test isgreen(dt[t₂])
    mark!(dt[t₂],2)
    @test isblue(dt[t₂])
    mark!(dt[t₂],3)
    @test isred(dt[t₂])

    # Mesh creation and basic properties. 
    vert = [0. 0.;1. 0.;1. 1.;0. 1.]'
    m0 = hpmesh(vert,sqrt(2)/2)
    @test length(m0.points) == 5
    @test length(m0.edgelist) == 8
    @test length(m0.trilist) == 4
    @test floattype(m0) == Float64
    @test inttype(m0) == Int32
    @test degtype(m0) == UInt8
    ∂m0 = BoundaryHPMesh(m0,1)
    @test ∂m0 == BoundaryHPMesh(m0,:dirichlet)
    @test length(edges(∂m0)) == 4
    vert1 = Float32.(vert)
    m1 = hpmesh(vert1,sqrt(2)/2)
    @test length(m1.points) == 5
    @test length(m1.edgelist) == 8
    @test length(m1.trilist) == 4
    @test floattype(m1) == Float32
    @test inttype(m1) == Int32
    @test degtype(m1) == UInt8
    setneumann!(x->x[2]==1,m1)
    ∂Dm1 = dirichletboundary(m1)
    ∂Nm1 = neumannboundary(m1)
    @test length(edges(∂Dm1)) == 3
    @test length(edges(∂Nm1)) == 1

    # Building mesh from points
    pts = [0. 0.;1. 0.;1. 1.;0. 1.]'
    T₁ = triangle(Int32[1,2,3],pts)
    T₂ = triangle(Int32[2,3,4],pts)
    @test isequal(longestedge(T₁), Edge{Int32}(1,3))
    eds₂ = tuple(edges(T₂)...)
    @test isequal(eds₂[1],Edge{Int32}(2,4))
    @test isequal(eds₂[2],Edge{Int32}(2,3))
    @test isequal(eds₂[3],Edge{Int32}(3,4))

    pts₂ = SVector{2,Float64}.([(0,0),(1,0),(0,1),(1,1)])
    T₃ = triangle(Int32[1,2,3],pts₂)
    T₄ = triangle(Int32[2,3,4],pts₂)
    @test isequal(T₁,T₃)
    @test isequal(T₂,T₄)
    trilist = Dictionary([T₃,T₄],[TriangleAttributes(),TriangleAttributes()])
    edgelist = Dictionary{Edge{Int32},EdgeAttributes{UInt8}}()
    for ed in edges(T₃)
        set!(edgelist,ed,EdgeAttributes{UInt8}(1,1,false))
    end
    for ed in edges(T₄)
        set!(edgelist,ed,EdgeAttributes{UInt8}(1,1,false))
    end
    mesh = HPMesh(pts₂,trilist,edgelist)
    @test degrees_of_freedom!(mesh)==length(pts₂)
    @test length(tagged_dof(mesh,1))==4
    @test !isempty(dof(mesh))
    empty!(dof(mesh))
    @test isempty(dof(mesh))

    @test edges(mesh) === keys(mesh.edgelist)

    edsmesh = mesh.edgelist
    setdegree!(edsmesh[longestedge(T₃)],4)
    @test !check_p_conformity(mesh)
    p_conformity!(mesh)
    @test check_p_conformity(mesh)
        
    @test repr(T₃) == "(2, 3, 1)"
    @test repr(TriangleAttributes()) == ":noref"

    @test repr(circmesh_graded_center(0.1,0.45)) == "HPMesh{Float64, Int32, UInt8}
  + 692 nodes.
  + 1320 triangles.
  + 2011 edges."
   
