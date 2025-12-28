using TrihpFEM
using Test
using Aqua

@testset "TrihpFEM.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(TrihpFEM)
    end
    @testset "Meshes" begin include("mesh.jl") end
end

