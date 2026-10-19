"""
   ref_integrate(p::ImmutablePolynomial)

Integrates `p` in the interval [-1,1]. The integration is performed exactly, with no quadratures.
"""
function ref_integrate(p::ImmutablePolynomial{F,X,N}) where {F,X,N}
    ip = Polynomials.integrate(p)
    ip(one(F))-ip(-one(F))
end
"""

    ref_integrate(p::PolyScalarField)

Integrates `p` in the reference triangle. The integration is performed exactly, with no quadratures.
"""
function ref_integrate(p::ProductPoly{2,F}) where F
    px,py = p.polys
    qy = Polynomials.integrate(py)
    x = ImmutablePolynomial((zero(F),one(F)),:x)
    qx = px*(qy(x)-qy(-one(F)))
    q = Polynomials.integrate(qx)
    q(one(F))-q(-one(F))
end
ref_integrate(p::PolySum) = ref_integrate(p.left) + ref_integrate(p.right)



"""

    ref_integrate(fun,sch)

Integrates `fun` in the reference triangle using the `Quadrature` `sch.  
"""
function ref_integrate(fun,sch::Quadrature{D,R,V,P}) where {D,R,V,P}
    (;weights,points)= sch
     #measure of the element. Weights for 1D quadratures are already scaled to [-1,1], so we multiply by 1. The reference triangle has area 2. 
    D*sum(w*fun(p) for (w,p) in zip(weights,points))
end


function edge_integrate₁(p::ProductPoly{2,F}) where F
    uno = one(F)
    px,py = p.polys
    i = py(-uno)
    i == zero(F) && return zero(F)
    qx = Polynomials.integrate(px)
    return i*(qx(uno)-qx(-uno))
end

function edge_integrate₂(p::ProductPoly{2,F}) where F
    uno = one(F)
    px,py = p.polys
    i = px(uno)
    i == zero(F) && return zero(F)
    qy = Polynomials.integrate(py)
    return i*(qy(uno)-qy(-uno))
end
function edge_integrate₃(p::ProductPoly{2,F}) where F
    px,py = p.polys
    qx = ImmutablePolynomial(py.coeffs)
    qq = Polynomials.integrate(px*qx)
    return sqrt(2)*(qq(one(F))-qq(-one(F)))
end

edge_integrate(p::ProductPoly{2},::Val{1}) = edge_integrate₁(p) 
edge_integrate(p::ProductPoly{2},::Val{2}) = edge_integrate₂(p) 
edge_integrate(p::ProductPoly{2},::Val{3}) = edge_integrate₃(p)

edge_integrate(p::PolySum,::Val{E}) where E = edge_integrate(p.left,Val(E)) + edge_integrate(p.right,Val(E))

