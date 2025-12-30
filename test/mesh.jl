    using ..TrihpFEM.Meshes
    using Test
    using Dictionaries

    #Edge creation and comparison.
    e₁ = Edge(3,4)
    e₂ = Edge((4,3))
    @test isequal(e₁,e₂)

    ea₁ = EdgeAttributes{UInt8}(0,1,false)

    de = Dictionary([e₁],[ea₁])
    @test e₁ in keys(de)
    @test de[e₂] == ae₁
    @test degree(de[e₁]) == 1
    @test !ismarked(de[e₁])
    @test marker(de[e₁]) == 1
    mark!(de[e₁])
    @test ismarked(de[e₂])
    setdegree!(de[e₁],3)
    @test typeof(de[e₁]) == UInt8
    @test degree(de[e₁]) == 3
    @test !isinterior(de[e₂])
    @test data(e₁) == (3,4)

    #Triangle creation and comparison
    t₁ = Triangle(1,2,3)
    t₂ = Triangle((2,3,1))
    @test isequal(e₁,e₂)

    t₃ = Triangle{UInt32}(1,3,2)
    @test isequal(t₁,t₃)

    ta₁ = TriangleAtttributes{UInt8,Float64}()
    ta₂ = TriangleAtttributes()
    @test ta₁==ta₂

    de = Dictionary([e₁],[ta₁])
    @test e₁ in keys(de)
    @test de[e₂] == ae₂
    @test degree(de[e₁]) == 1
    @test !ismarked(de[e₂])
    @test marker(de[e₁]) == 1
    mark!(de[e₁])
    @test ismarked(de[e₂])
    setdegree!(de[e₁],3)
    @test typeof(de[e₁]) == UInt8
    @test degree(de[e₁]) == 3
    @test !isinterior(de[e₂])
    @test data(e₁) == (3,4)

    # Mesh creation and basic properties. 
    vert = [0. 0.;1. 0.;1. 1.;0. 1.]'
    m0 = hpmesh(vert,sqrt(2)/2)
    @test length(m0.points) == 5
    @test length(m0.edgelist) == 8
    @test length(m0.trilist) == 4
    @test floattype(m0) == Float64
    @test inttype(m0) == Int32
    @test degtype(m0) == UInt8
    ∂m0 = BoundaryHPEdge(m0,1)
    @test ∂m0 == BoundaryHPEdge(m0,:dirichlet)
    @test length(edges(∂m0)) == 4
    vert1 = Float32.(vert)
    m1 = hpmesh(vert,sqrt(2)/2)
    @test length(m1.points) == 5
    @test length(m1.edgelist) == 8
    @test length(m1.trilist) == 4
    @test floattype(m1) == Float32
    @test inttype(m1) == Int32
    @test degtype(m1) == UInt8

    
    
    

    
    
    
    
