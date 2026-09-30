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



