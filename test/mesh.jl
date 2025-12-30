    using ..TrihpFEM.Meshes
    using Test
    using Dictionaries

    #Edge creation and comparison.
    e₁ = Edge(3,4)
    e₂ = Edge((4,3))
    @test isequal(e₁,e₂)
    ea₁ = Meshes.EdgeAttributes{UInt8}(0,1,false)

    de = Dictionary([e₁],[ea₁])
    @test e₁ in keys(de)
    @test de[e₂] == ae₁
    @test degree(de[e₁]) == 1
    @test Meshes.!ismarked(de[e₁])
    @test Meshes.marker(de[e₁]) == 1
    Meshes.mark!(de[e₁])
    @test Meshes.ismarked(de[e₂])
    Meshes.setdegree!(de[e₁],3)
    @test typeof(de[e₁].marker[]) == UInt8
    @test degree(de[e₁]) == 3
    @test Meshes.!isinterior(de[e₂])
    @test meshes.data(e₁) == (3,4)

    #Triangle creation and comparison
    t₁ = Triangle(1,2,3)
    t₂ = Triangle((2,3,1))
    @test isequal(e₁,e₂)

    t₃ = Triangle{UInt32}(1,3,2)
    @test isequal(t₁,t₃)

    ta₁ = Meshes.TriangleAtttributes{UInt8,Float64}()
    ta₂ = Meshes.TriangleAtttributes()
    @test ta₁==ta₂

    dt = Dictionary([t₁],[ta₁])
    @test t₂ in keys(dt)
    @test Meshes.!ismarked(dt[e₂])
    Meshes.mark!(dt[e₁],1)
    @test Meshes.ismarked(dt[e₂])
    @test Meshes.isgreen(dt[e₂])
    Meshes.mark!(dt[e₂],2)
    @test Meshes.isblue(dt[e₂])
    Meshes.mark!(dt[e₂],3)
    @test Meshes.isred(dt[e₂])

    # Mesh creation and basic properties. 
    vert = [0. 0.;1. 0.;1. 1.;0. 1.]'
    m0 = hpmesh(vert,sqrt(2)/2)
    @test length(m0.points) == 5
    @test length(m0.edgelist) == 8
    @test length(m0.trilist) == 4
    @test Meshes.floattype(m0) == Float64
    @test Meshes.inttype(m0) == Int32
    @test Meshes.degtype(m0) == UInt8
    ∂m0 = BoundaryHPEdge(m0,1)
    @test ∂m0 == BoundaryHPEdge(m0,:dirichlet)
    @test length(edges(∂m0)) == 4
    vert1 = Float32.(vert)
    m1 = hpmesh(vert,sqrt(2)/2)
    @test length(m1.points) == 5
    @test length(m1.edgelist) == 8
    @test length(m1.trilist) == 4
    @test Meshes.floattype(m1) == Float32
    @test Meshes.inttype(m1) == Int32
    @test Meshes.degtype(m1) == UInt8

    
    
    

    
    
    
    
