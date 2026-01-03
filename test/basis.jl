using ..TrihpFEM.Poly
using Test
using StaticArrays

@test_throws StandardBasis((1,2,4)) ArgumentError("Degrees does not satisfy p conformity.")

B = StandardBasis((2,2,3))
C = StandardBasis((2,2,3))
res0 = [1.0  0.0  -0.5  0.0;
           0.0  0.0  -0.0  0.0;
           -0.5 -0.0 0.25  -0.0;
           0.0  0.0  -0.0  0.0]
res023077 = [1.0         0.23        -0.42065     -0.314583;
              0.77        0.1771      -0.3239      -0.242229;
              0.38935     0.0895505   -0.16378     -0.122483;
             -0.0136675  -0.00314352   0.00574923   0.00429956]

k = 0
for (i,b) in enumerate(B)
    @test res0[i] ≈ b([0.,0.])
    @test res023077[i] ≈ b(SVector([0.23,0.77]))
    k = k+1
end
@test length(B) == k