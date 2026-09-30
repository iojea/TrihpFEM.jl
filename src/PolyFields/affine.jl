
###########################################################################################
#########################           AffineToRef        ###########################
###########################################################################################
"""

    AffineToRef{F}

A struct for defining and updating an affine transformation from the reference triangle to some other triangle. 
"""
struct AffineToRef{D,F <: Number}
    A::Tensor{2,D,F}
    b::Tensor{1,D,F}
    function AffineToRef(vert)
        D = length(vert)-1
        D ∈ (1,2) || throw(ArgumentError("Two or three vertices are needed."))
        A = affinetoref_matrix(Val(D),vert)
        b = affinetoref_vec(Val(D),vert)
        F = promote_type(eltype(A),eltype(b))
        A = Tensor{2,D}((F(x) for x in A)...)
        b = Tensor{1,D}((F(x) for x in b)...)
        return new{D,F}(A,b)
    end
end
function affinetoref_matrix(::Val{2}, vert)
    return Tensor{2, 2}(t -> (vert[t[2]+1][t[1]] - vert[t[2]][t[1]])/2)
end
# function affinetoref_matrix(::Val{1}, vert)
#     return Tensor{2,1}(i -> (last(vert)[i] - first(vert)[i])/2)
# end

function affinetoref_vec(::Val{D},vert) where D
    return Tensor{1, D}(i -> (first(vert)[i] + last(vert)[i]) / 2)
end

(aff::AffineToRef)(x)  = aff.A * x + aff.b


jac(aff::AffineToRef{2,F}) where F = abs(det(aff.A))
jac(aff::AffineToRef{1,F}) where F = norm(aff.A)
# area(x,y,z) = 0.5abs(x[1]*(y[2]-z[2])+y[1]*(z[2]-x[2])+z[1]*(x[2]-y[2]))
# area(v::Vector) = area(v...)
area(t::AffineToRef{2,F}) where F = 2jac(t)


# A trait for evaluation of Field
abstract type EvalType end
struct Eval <: EvalType end
struct Compose <: EvalType end
struct Pass <: EvalType end

evaltype(_) = Pass()
evaltype(::Function) = Compose()
evaltype(::PolyScalarField) = Eval()

evaluate(::Eval, f, ::AffineToRef, x) = f(x)
evaluate(::Compose, f, t::AffineToRef, x) = f(t(x))
evaluate(::Pass, f, ::AffineToRef, _) = f
