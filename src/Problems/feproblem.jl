"""

    FEProblem(a,b,space,g)

Defines a Finite Element Problem (no refinement) of the form `a(u,v)=b(v)` over the space `space` with Dirichlet boundary data `g`.

# Arguments
+ `a`: bilinear form.
+ `b`: linear form.
+ `g` a function defining the Dirichlet boundary condition.
"""

_toform(a::Form) = a
_toform(a::Term{C,O,T,N,M}) where {C,O,T,N,M} = Form{N}((a,))
struct FEProblem{F₁,F₂}
    a::F₁
    b::F₂
    g
    K
    rhs
    function FEProblem(a,b,g)
        K = assembly_matrix(a)
        rhs = assembly_rhs(b)
        fa = _toform(a)
        fb = _toform(b)
        new{typeof(fa),typeof(fb)}(fa,fb,g,K,rhs)
    end        
end
FEProblem(a,b) = FEProblem(a,b,x->zero(eltype(x)))

function solve(prob::FEProblem{F₁,F₂}) where {F₁,F₂}
    mesh = domainmesh(first(prob.a.terms).measure)
   (;points,edgelist,dofs) = mesh
    vals = FixedSizeVector{floattype(mesh)}(undef,ndof(mesh))
    dirichlet = Vector{inttype(mesh)}(undef,0)
    for e in keys(edgelist)
        if istagged(edgelist[e],1)
            locdof = dofs.by_edge[e]
            i1 = first(locdof)
            i2 = last(locdof)
            x1 = points[i1]
            x2 = points[i2]
            n = length(locdof)
            for (i,k) in enumerate(locdof)
                vals[k] = prob.g(((n-i)*x1+(i-1)*x1)/(n-1))
            end
            push!(dirichlet,locdof...)
        end
    end
    unique!(dirichlet)
    valid = setdiff(1:ndof(mesh),dirichlet)
    vals[valid] = prob.K[valid,valid]\prob.rhs[valid]
    return FESolution(mesh,vals)
end 

# function FEProblem(a,b,space,g)
#     matrix = integrate(a,space)
#     rhs = integrate(b,space)
# end

