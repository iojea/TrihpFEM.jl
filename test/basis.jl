using ..TrihpFEM.Poly
using Test
using StaticArrays

@test_throws ArgumentError("Degrees does not satisfy p conformity.") StandardBasis((1,2,4))

B = StandardBasis((2,3,4))
res0 =  [1.0   0.0  -0.5    0.0;
         0.0   0.0  -0.0    0.0;
        -0.5  -0.0   0.25  -0.0]

res023077 = [1.0       0.77     0.38935    -0.0136675;
             0.23      0.1771   0.0895505  -0.00314352;
            -0.42065  -0.3239  -0.16378     0.00574923]
k = 0
for b in B
    (i,j) = degs(b)
    @test res0[i+1,j+1] ≈ b([0.,0.])
    @test res023077[i+1,j+1] ≈ b(SVector(0.23,0.77))
    global k = k+1
end

@test length(B) == k