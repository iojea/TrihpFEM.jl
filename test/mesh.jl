    using ..TrihpFEM.Meshes
    using Test
    using Dictionaries
    using StaticArrays

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
    @test typeof(Meshes.tag(de[e₁])) == UInt8
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
    ∂m0 = BoundaryHPMesh(m0,1)
    @test ∂m0 == BoundaryHPMesh(m0,:dirichlet)
    @test length(edges(∂m0)) == 4
    vert1 = Float32.(vert)
    m1 = hpmesh(vert1,sqrt(2)/2)
    @test length(m1.points) == 5
    @test length(m1.edgelist) == 8
    @test length(m1.trilist) == 4
    @test Meshes.floattype(m1) == Float32
    @test Meshes.inttype(m1) == Int32
    @test Meshes.degtype(m1) == UInt8
    Meshes.set_neumann!(x->x[2]==1,m1)
    ∂Dm1 = dirichletboundary(m1)
    ∂Nm1 = neumannboundary(m1)
    @test length(edges(∂Dm1)) == 3
    @test length(edges(∂Nm1)) == 1

    # Building mesh from points
    pts = [0. 0.;1. 0.;1. 1.;0. 1.]'
    T₁ = triangle(Int32[1,2,3],pts)
    T₂ = triangle(Int32[2,3,4],pts)
    @test isequal(Meshes.longestedge(T₁), Edge{Int32}(1,3))
    eds₂ = tuple(edges(T₂)...)
    @test isequal(eds₂[1],Edge{Int32}(2,4))
    @test isequal(eds₂[2],Edge{Int32}(2,3))
    @test isequal(eds₂[3],Edge{Int32}(3,4))
        

    
    
    

    
    
    
    
