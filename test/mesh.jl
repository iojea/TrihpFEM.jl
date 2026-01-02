    using ..TrihpFEM.Meshes
    using Test
    using Dictionaries
    using StaticArrays
    using Makie

    #Edge creation and comparison.
    e₁ = Edge{Int32}(UInt8(3),4)
    e₂ = Edge{Int32}((4,3))
    e₃ = Edge{Int32}(SVector(4,3))
    e₄ = Edge{Int32}(3,4)
    e₅ = Edge{Int32}(e for e in e₄)
    e₆ = Edge{Int32}([4,3])
    @test isequal(e₁,e₂) && isequal(e₁,e₃) && isequal(e₄,e₅)
    @test isequal(e₁,e₆) && isequal(e₄,e₆)
    @test typeof(e₄) == Edge{Int32}
    @test isequal(e₁,Edge(3,4))
    @test isequal(e₂,Edge((3,4)))
    @test isequal(e₃,Edge(SVector(4,3)))
    @test isequal(e₄,Edge(3,4))
    @test isequal(e₅,Edge(e for e in e₅))
    @test isequal(e₆,Edge([4,3]))
    
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
    @test !isboundary(de[e₂])
    @test data(e₁) == (3,4)

    #Triangle creation and comparison
    t₁ = Triangle{Int32}(Int8(1),2,3)
    t₂ = Triangle{Int32}((2,3,1))
    t₃ = Triangle{Int32}(SVector(3,1,2))
    t₄ = Triangle{UInt32}(2,1,3)
    t₅ = Triangle{Int32}(t for t in t₄)
    t₆ = Triangle{Int32}([1,2,3])
    @test isequal(t₁,t₂) && isequal(t₁,t₃) && isequal(t₁,t₄)
    @test isequal(t₄,t₅) && isequal(t₄,t₆)
    @test eltype(data(t₄)) == UInt32
    @test isequal(t₁,Triangle(1,2,3))
    @test isequal(t₂,Triangle((2,3,1)))
    @test isequal(t₃,Triangle(SVector(3,1,2)))
    @test isequal(t₄,Triangle(2,1,3))
    @test isequal(t₅,Triangle(t for t in t₅))
    @test isequal(t₆,Triangle([1,2,3]))

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
    @test repr(m0.edgelist) == "EdgeList{Int32,UInt8} with 8 elements"
    @test repr(m0.trilist) == "TriangleList{Int32,UInt8,Float64} with 4 elements"
    @test length(m0.points) == 5
    @test length(m0.edgelist) == 8
    @test length(m0.trilist) == 4
    @test floattype(m0) == Float64
    @test inttype(m0) == Int32
    @test degtype(m0) == UInt8
    ∂m0 = BoundaryHPMesh(m0,1)
    @test ∂m0 == BoundaryHPMesh(m0,:dirichlet)
    @test length(edges(∂m0)) == 4
    @test domainmesh(∂m0) === m0
    @test domainmesh(domainmesh(∂m0)) === m0
    vert1 = Float32.(vert)
    m1 = hpmesh(vert1,sqrt(2)/2)
    @test length(m1.points) == 5
    @test length(m1.edgelist) == 8
    @test length(m1.trilist) == 4
    @test floattype(m1) == Float32
    @test inttype(m1) == Int32
    @test degtype(m1) == UInt8
    setneumann!(x->x[1]*x[2]==0,m1)
    ∂Dm1 = dirichletboundary(m1)
    ∂Nm1 = neumannboundary(m1)
    @test length(edges(∂Dm1)) == 2
    @test length(edges(∂Nm1)) == 2
    setboundary!(:dirichlet,x->x[2]==0,m1)
    ∂Dm1 = dirichletboundary(m1)
    ∂Nm1 = neumannboundary(m1)
    @test length(edges(∂Dm1)) == 3
    @test length(edges(∂Nm1)) == 1
    setdirichlet!(x->x[1]==0,m1)
    ∂Dm1 = dirichletboundary(m1)
    ∂Nm1 = neumannboundary(m1)
    @test length(edges(∂Dm1)) == 4
    @test length(edges(∂Nm1)) == 0
    
    

    # Building mesh from points
    pts = Float32[0. 0.;1. 0.;1. 1.;0. 1.]'
    T₁ = triangle(Int16[1,2,3],pts)
    T₂ = triangle(Int16[2,3,4],pts)
    @test isequal(longestedge(T₁), Edge{Int16}(1,3))
    eds₂ = tuple(edges(T₂)...)
    @test isequal(eds₂[1],Edge{Int16}(2,4))
    @test isequal(eds₂[2],Edge{Int16}(2,3))
    @test isequal(eds₂[3],Edge{Int16}(3,4))
    mmtris = Dictionary([T₁,T₂],[TriangleAttributes{UInt8,Float32}(0,0,0) for _ in 1:2])
    mmedges = Dictionary{Edge{Int16},EdgeAttributes{UInt8}}()
    for ed in edges(T₁)
        set!(mmedges,ed,EdgeAttributes{UInt8}(1,1,false))
    end
    for ed in edges(T₂)
        set!(mmedges,ed,EdgeAttributes{UInt8}(1,1,false))
    end
    mm = HPMesh(pts,mmtris,mmedges)
    @test typeof(mm) == HPMesh{Float32,Int16,UInt8}

    pts₂ = SVector{2,Float64}.([(0,0),(1,0),(0,1),(1,1)])
    T₃ = triangle(Int32[1,2,3],pts₂)
    T₄ = triangle(Int32[2,3,4],pts₂)
    @test isequal(T₁,T₃)
    @test isequal(T₂,T₄)
    trilist = Dictionary([T₃,T₄],[TriangleAttributes(),TriangleAttributes()])
    edgelist = Dictionary{Edge{Int32},EdgeAttributes{UInt8}}()
    @test repr(edgelist) == "EdgeList{Int32,UInt8} with 0 elements"
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
        
    @test repr(TriangleAttributes()) == ":noref"

    # REFINE h
    # This test should be replaced by something else that tests h-refinement. For the moment, it is a quick way to test a large number of functions. 
    @test typeof(circmesh_graded_center(0.1,0.45)) == HPMesh{Float64,Int32,UInt8}
    @test typeof(squaremesh(1,0.1)) == HPMesh{Float64,Int32,UInt8}
    cm = circmesh(1,0.1)

    #testing plots
    # marking for refinement. 
    for t in pairs(cm.trilist)
        if rand()<0.1
        mark!(last(t))
        mark!.(getindices(cm.edgelist,edges(first(t))))
        end
    end
    Meshes._h_conformity!(cm)
    plt = plothpmesh(cm)
    @test plt isa Makie.FigureAxisPlot

    for e in cm.edgelist
       setdegree!(e,rand(1:15))
    end
    plt2 = degplot(cm)
    @test plt isa Makie.FigureAxisPlot

    @test typeof(cm) == HPMesh{Float64,Int32,UInt8}
    @test contains(repr(cm),"HPMesh{Float64, Int32, UInt8}")
    long_print = sprint((io, x) -> show(IOContext(io, :limit => true), MIME("text/plain"), x), cm)
    @test contains(long_print,"HPMesh{Float64, Int32, UInt8}")
    @test contains(long_print,"Dictionary{Triangle{Int32}, TriangleAttributes{UInt8, Float64}")
    @test contains(long_print,"Dictionary{Edge{Int32}, EdgeAttributes{UInt8}}")

    pts₃ = [0. 0.;1. 0.;1. 0.5;1. 1.;0. 1.;0. 0.5;0.25 0.25;0.5 0.25;0.5 0.5;0.25 0.5]';
    segs₃ = [1 2;2 3;3 4;4 5;5 6;6 1;7 8;8 9;9 10;10 7]';
    tags₃ = [:dirichlet,:dirichlet,:neumann,:neumann,:neumann,
             :dirichlet,:dirichlet,:dirichlet,:neumann,:neumann];
    hole = [0.3 0.3]';
    mhole = hpmesh(pts₃,0.1;segments=segs₃,tags=tags₃,holes=hole)
    @test typeof(mhole) == HPMesh{Float64,Int32,UInt8}