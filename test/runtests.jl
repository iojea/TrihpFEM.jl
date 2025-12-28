using TrihpFEM
using Test
using Aqua

@testset "TrihpFEM.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(TrihpFEM)
    end
    # Write your tests here.
end
