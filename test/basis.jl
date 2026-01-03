using ..TrihpFEM.Poly
using Test
using StaticArrays

@test_throws StandardBasis((1,2,4)) ArgumentError("Degrees does not satisfy p conformity.")
B = StandardBasis((2,2,3))
C = StandardBasis((2,2,3))
results = [1.0  0.0  -0.5  0.0;
           0.0  0.0  -0.0  0.0;
           -0.5 -0.0 0.25  -0.0;
           0.0  0.0  -0.0  0.0]
for (i,b) in enumerate(B)
    for(j,c) in enumerate(B)
        @test b([0,0]