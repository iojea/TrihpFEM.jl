
###########################################################################################
#########################           AffineToRef        ###########################
###########################################################################################
"""

    AffineToRef{F}

A struct for defining and updating an affine transformation from the reference triangle to some other triangle. 
"""
struct AffineToRef{D,F <: Real,L}
    A::Tensor{2,D,F,L,Tuple{D,D}}
    b::Tensor{1,D,F,D,Tuple{D}}
    function AffineToRef(vert)
        D = length(first(vert)) 
        length(vert)==D+1 || throw(ArgumentError("Two or three vertices are needed."))
        A = affinetoref_matrix(Val(D),vert)
        b = affinetoref_vec(Val(D),vert)
        F = promote_type(eltype(A),eltype(b))
        A = Tensor{2,D}((F(x) for x in A)...)
        b = Tensor{1,D}((F(x) for x in b)...)
        L = D*D
        return new{D,F,L}(A,b)
    end
end
function affinetoref_matrix(::Val{2},v)
    v1,v2,v3 = v
    tup = ((v2[1]-v1[1])/2,(v2[2]-v1[2])/2,(v3[1]-v2[1])/2,(v3[2]-v2[2])/2)
    Tensor{2,2,eltype(tup),4,Tuple{2,2}}(tup)
end
# function affinetoref_matrix(::Val{1}, vert)
#     return Tensor{2,1}(i -> (last(vert)[i] - first(vert)[i])/2)
# end

function affinetoref_vec(::Val{D},v) where D
    @inbounds tup = ntuple(i->(first(v)[i]+last(v)[i])/2,D)
    T = eltype(tup)
    Tensor{1,D,T,D,Tuple{D}}(tup)
end


AffineToRef(tri::Triangle,mesh::HPMesh) = AffineToRef(mesh.points[tri])

(aff::AffineToRef)(x)  = aff.A * x + aff.b


jac(aff::AffineToRef{2,F}) where F = abs(det(aff.A))
jac(aff::AffineToRef{1,F}) where F = norm(aff.A)
# area(x,y,z) = 0.5abs(x[1]*(y[2]-z[2])+y[1]*(z[2]-x[2])+z[1]*(x[2]-y[2]))
# area(v::Vector) = area(v...)
area(t::AffineToRef{2,F}) where F = 2jac(t)


# A trait for evaluation of Field
# abstract type EvalType end
# struct Eval <: EvalType end
# struct Compose <: EvalType end
# struct Pass <: EvalType end

# evaltype(_) = Pass()
# evaltype(::Function) = Compose()
# evaltype(::PolyScalarField) = Eval()

# evaluate(::Eval, f, ::AffineToRef, x) = f(x)
# evaluate(::Compose, f, t::AffineToRef, x) = f(t(x))
# evaluate(::Pass, f, ::AffineToRef, _) = f

# @inline function LinearAlgebra.inv(t::Tensor{2,2})
#     id = 1/det(t)
#     @inbounds tup = promote(t[4]*id,-t[2]*id,-t[3]*id,t[1]*id)
#     Tensor{2,2,eltype(tup),4,Tuple{2,2}}(tup)
# end
