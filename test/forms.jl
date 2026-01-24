using Test

A = rand(2,2)
f(x) = x[1]*x[2]
@form a(u,v) = ∫((A*∇(u))⋅∇(v))*dΩ
@form b(v) = ∫(f*v)*dΩ

@test typeof(a)<:Form{2}
@test typeof(b)<:Form{1}


