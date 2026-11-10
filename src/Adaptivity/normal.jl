struct Normal end

function _normal(p::Point2D,q::Point2D)
    r = q-p
    R = norm(r)
    return Point2D(r[2]/R,-r[1]/R)
end

function (::Normal)(e::Edge,bm::BoundaryMeasure)
    tag(e)==bm.kind || throw(ArgumentError("Edge $e does not belong to the measure $bm."))
    msh = mesh(bm)
    adj = adjacents(e)
    tri = Triangle(e[1],e[2],adj[1])
    p1,p2,p3 = view(msh.points,tri)
    if orient(p1,p2,p3)
        return _normal(p1,p2)
    else
        return _normal(p2,p1)
    end
end

# function (::Normal)(tri::Triangle,e::Edge,bm::SkeletonMeasure)
#     nothing ## Check whether is better to pass the edge as a parameter or just an index.
# end
