module MeshTests
    using TrihpFEM.Meshes
    using Test
    using Dictionaries

    p = Point2(0.,1)
    q = Point2{Float64}(0,1)
    r = Point2((0.,1))
    @test p==q
    @test p==r
    @test q==r
    @test zero(p)==Point2(0.,0.)

    pf = Point2f(1,0.)
    qf = Point2f((1.,0))
    rf = Point2d(pf)
    @test pf==qf
    @test pf==rf
    @test qf==rf
    @test zero(pf)==Point2f(0,0)

    e₁ = Edge(3,4)
    e₂ = Edge((4,3))
    @test isequal(e₁,e₂)

    ea₁ = EdgeAttributes{UInt8}(0,1,false)
    ea₂ = EdgeAttributes(UInt8(0),UInt8(1),false)
    @test ea₁==ea₂

    t₁ = Triangle(1,2,3)
    t₂ = Triangle((2,3,1))
    @test isequal(e₁,e₂)

    t₃ = Triangle{UInt32}(1,3,2)
    @test isequal(t₁,t₃)

    de = Dictionary([e₁],[ea₁])
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

end;