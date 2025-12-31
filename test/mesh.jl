    using ..TrihpFEM.Meshes
    using Test
    using Dictionaries

    #Edge creation and comparison.
    e₁ = Edge(3,4)
    e₂ = Edge((4,3))
    @test isequal(e₁,e₂)
    ea₁ = Meshes.EdgeAttributes{UInt8}(1,0,false)

    de = Dictionary([e₁],[ea₁])
    @test e₁ in keys(de)
    @test de[e₂] == ea₁
    @test degree(de[e₁]) == 1
    @test !Meshes.ismarked(de[e₁])
    @test Meshes.tag(de[e₁]) == 0
    Meshes.mark!(de[e₁])
    @test Meshes.ismarked(de[e₂])
    Meshes.setdegree!(de[e₁],3)
    @test typeof(tag(de[e₁])) == UInt8
    @test degree(de[e₁]) == 3
    @test Meshes.isinterior(de[e₂])
    @test Meshes.data(e₁) == (3,4)

    #Triangle creation and comparison
    t₁ = Triangle(1,2,3)
    t₂ = Triangle((2,3,1))
    @test isequal(e₁,e₂)

    t₃ = Triangle{UInt32}(1,3,2)
    @test isequal(t₁,t₃)

    ta₁ = Meshes.TriangleAttributes{UInt8,Float64}()
    ta₂ = Meshes.TriangleAttributes()
    @test all(getproperty(ta₁,prop)[]==getproperty(ta₂,prop)[] for prop in propertynames(ta₁))

    dt = Dictionary([t₁],[ta₁])
    @test t₂ in keys(dt)
    @test !Meshes.ismarked(dt[t₂])
    Meshes.mark!(dt[t₁],1)
    @test Meshes.ismarked(dt[t₂])
    @test Meshes.isgreen(dt[t₂])
    Meshes.mark!(dt[t₂],2)
    @test Meshes.isblue(dt[t₂])
    Meshes.mark!(dt[t₂],3)
    @test Meshes.isred(dt[t₂])

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

    
    
    

    
    
    
    
