"""

    ref_integrate(p::PolyScalarField)

Integrates `p` in the reference triangle. The integration is performed exactly, with no quadratures.
"""
function ref_integrate(p::BiPoly{F,X,Y}) where {F,X,Y}
    (;px,py) = p
    qy = Polynomials.integrate(py)
    x = ImmutablePolynomial((zero(F),one(F)),X)
    qx = px*(qy(x)-qy(-one(F)))
    q = Polynomials.integrate(qx)
    q(one(F))-q(-one(F))
end
ref_integrate(p::PolySum) = ref_integrate(p.left) + ref_integrate(p.right)



"""

    ref_integrate(fun,sch)

Integrates `fun` in the reference triangle using the `Quadrature` `sch.  
"""
function ref_integrate(fun,sch)
    (;weights,points)= sch
    2sum(w*fun(p) for (w,p) in zip(weights,points))
end

