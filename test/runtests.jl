using TrihpFEM
using Test
using Aqua

@testset "TrihpFEM.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(TrihpFEM,stale_deps=false)
    end
    @testset "Meshes" begin include("mesh.jl") end
    @testset "Polys" begin include("poly.jl") end
    @testset "Basis" begin include("basis.jl") end
end

