affine(tri::Triangle, m::DomainMeasure) = AffineToRef(mesh(m).points[tri])

function affine(e::Edge, m::BoundaryMeasure)
    eds = elements(m)
    adj = adjacents(eds[e])
    adj[2] == 0 || throw(ArgumentError("Is not possible to create an `AffineToRef` transformation from an interior edge."))
    tri = Triangle(e[1], e[2], adj[1])
    return affine(tri, m)
end
