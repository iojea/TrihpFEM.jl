"""
$(SIGNATURES)

Is a cache for the refining process. It stores data that will be updated for each triangle.
"""
struct RefineAux{I <: Integer, P <: Integer}
    i::Base.RefValue{I}
    degs::MVector{6, P}
    dots::MVector{6, I}
end

function RefineAux{I, P}() where {I,P}
    return RefineAux(MVector{6, P}(zeros(6)), MVector{6, I}(zeros(6)))
end

function RefineAux(i, mesh::HPMesh{F, I, P}) where {F,P,I}
    return RefineAux(Ref(i), MVector{6, P}(zeros(6)), MVector{6, I}(zeros(6)))
end

"""
    $(SIGNATURES)

Marks triangles with vertices `vert` for which `estim(vert,estim_param)` returns `true`, and then run the `_h_conformity!` routine in order to propagate the markings to adjacent triangles as needed. 
"""
function mark!(estim::Function, mesh::HPMesh{F, I, P}; estim_params...) where {F, I, P}
    (; points, trilist, edgelist) = mesh
    tri = MMatrix{2, 3}(zero(F) for i in 1:2,j in 1:3)
    for t in triangles(trilist)
        for (i, tt) in enumerate(t)
            tri[:, i] .= points[tt]
        end
        if estim(tri; estim_params...)
            mark!.(getindices(edgelist, edges(t))) #AQuí había un FOR que cambié
        end
    end
    return _h_conformity!(mesh)
end

"""
    $(SIGNATURES)

checks if a point `p` belongs to the triangle with vertices `a`,`b` and `c`. 
"""
function intriangle(p::T, a::V, b::V, c::V) where {T <: AbstractArray, V <: AbstractArray}
    x = orient(a, b, p)
    y = orient(b, c, p)
    z = orient(c, a, p)
    if abs(x+y+z) == 3
        return 1
    elseif x*y*z == 0
        return 0
    else
        return -1
    end
end 
intriangle(p::T, vert::M) where {T <: AbstractArray, M <: AbstractArray} = intriangle(p, eachcol(vert)...)


"""
    $(SIGNATURES)

Performs the marking of previously un-marked triangles in order to avoid hanging nodes.
"""

function _h_conformity!(mesh::HPMesh)
    (; trilist, edgelist) = mesh
    still = true
    while still
        still = false
        for t in triangles(trilist)
            num_marked = count(ismarked, getindices(edgelist, edges(t)))
            if num_marked > 0
                long_edge = edgelist[longestedge(t)]
                if !ismarked(long_edge)
                    mark!(long_edge)
                    still = true
                    num_marked += 1
                end
                mark!(trilist[t], num_marked)
            end
        end
    end
    return
end


"""
  $(SIGNATURES)

performs the refinement of Red marked triangles.   
"""
function refine_red!(t::Triangle{I}, mesh::HPMesh{F, I, P}, refaux::RefineAux{I, P}) where {F,I,P}
    (; points, edgelist, trilist) = mesh
    (; i, degs, dots) = refaux
    dots[1:3] .= t
    t_edges = edges(t)
    degs[1:3] .= (degree(edgelist[e]) for e in t_edges)
    degs[4:6] .= max.(abs.(degs[SVector(1, 2, 3)] - degs[SVector(3, 1, 2)]), one(P))
    for (j,edge) in enumerate(t_edges)
        ep = edgelist[edge]
        middle = seen(ep)
        if middle == 0
            points[i[]] = SVector(sum(points[edge]) / 2)
            dots[j+3] = i[]
            seen!(ep,i[])
            m = tag(edgelist[edge])
            set!(edgelist, Edge(edge[1], i[]), EdgeAttributes(I,degs[j], m, false))
            set!(edgelist, Edge(i[], edge[2]), EdgeAttributes(I,degs[j], m, false))
            i[] += 1
        else
            dots[j+3] = middle
            seen!(ep,zero(I))
        end
    end
    for (j,edge) in enumerate(t_edges)
        s = isseen(edgelist[edge])
        setadjacent!(edgelist[Edge(edge[1],dots[j+3])],1+s,dots[3+mod1(j-1,3)])
        setadjacent!(edgelist[Edge(dots[j+3],edge[2])],1+s,dots[3+mod1(j-2,3)])
    end
    set!(edgelist, Edge(dots[SVector(6, 4)]), EdgeAttributes(I,degs[4], zero(P), false))
    set!(edgelist, Edge(dots[SVector(4, 5)]), EdgeAttributes(I,degs[5], zero(P), false))
    set!(edgelist, Edge(dots[SVector(5, 6)]), EdgeAttributes(I,degs[6], zero(P), false))
    set!(trilist, Triangle(dots[SVector(1, 4, 6)]), TriangleAttributes{F,P}())
    set!(trilist, Triangle(dots[SVector(4, 2, 5)]), TriangleAttributes{F,P}())
    set!(trilist, Triangle(dots[SVector(6, 5, 3)]), TriangleAttributes{F,P}())
    set!(trilist, Triangle(dots[SVector(5, 6, 4)]), TriangleAttributes{F,P}())
    nothing
end

"""
  $(SIGNATURES)

performs the refinement of Blue marked triangles.   
"""
function refine_blue!(t::Triangle{I}, mesh::HPMesh{F, I, P}, refaux::RefineAux{I, P}) where {F <: AbstractFloat, I <: Integer, P <: Integer}
    (; points, edgelist, trilist) = mesh
    (; i, degs, dots) = refaux
    dots[1:3] .= t
    t_edges = edges(t)
    degs[1:3] .= (degree(edgelist[e]) for e in t_edges)
    degs[4] = max(maximum(abs, degs[SVector(1, 3)] - degs[SVector(2, 1)]), one(P))
    if ismarked(edgelist[t_edges[2]])
        degs[5] = max(maximum(abs, degs[SVector(1, 2)] - degs[SVector(2, 4)]), one(P))
        for j in 1:2
            edge = t_edges[j]
            ep = edgelist[edge]
            middle = seen(ep)
            if middle > 0
                dots[j + 3] = middle
                seen!(ep,zero(I))
            else
                points[i[]] = SVector(sum(points[edge]) / 2)
                dots[j + 3] = i[]
                seen!(ep,i[])
                m = tag(edgelist[edge])
                set!(edgelist, Edge(edge[1], i[]), EdgeAttributes(I,degs[j], m, false))
                set!(edgelist, Edge(i[], edge[2]), EdgeAttributes(I,degs[j], m, false))
                i[] += 1
            end
        end
        for j in 1:2
            edge = t_edges[j]
            s = isseen(edgelist[edge])
            setadjacent!(edgelist[Edge(edge[1],dots[j+3])],1+s,dots[j+2])
            setadjacent!(edgelist[Edge(dots[j+3],edge[2])],1+s,dots[6-j])
        end
        set!(edgelist, Edge(dots[SVector(3, 4)]), EdgeAttributes(I,degs[4], zero(P), false))
        set!(edgelist, Edge(dots[SVector(5, 4)]), EdgeAttributes(I,degs[5], zero(P), false))
        set!(trilist, Triangle(dots[SVector(1, 4, 3)]), TriangleAttributes{F,P}())
        set!(trilist, Triangle(dots[SVector(4, 2, 5)]), TriangleAttributes{F,P}())
        set!(trilist, Triangle(dots[SVector(4, 5, 3)]), TriangleAttributes{F,P}())
    elseif ismarked(edgelist[t_edges[3]])
        degs[5] = max(maximum(abs, degs[SVector(1, 3)] - degs[SVector(3, 4)]), one(P))
        for j in 0:1
            edge = t_edges[1 + 2j]
            ep = edgelist[edge]
            middle = seen(ep)
            if middle > 0
                dots[j + 4] = middle
                seen!(ep,zero(I))
            else
                points[i[]] = SVector(sum(points[edge]) / 2)
                dots[j + 4] = i[]
                seen!(ep, i[])
                m = tag(edgelist[edge])
                set!(edgelist, Edge(edge[1], i[]), EdgeAttributes(I,degs[1 + 2j], m, false))
                set!(edgelist, Edge(i[], edge[2]), EdgeAttributes(I,degs[1 + 2j], m, false))
                i[] += 1
            end
        end
        for j in 0:1
            edge = t_edges[1+2j]
            s = isseen(edgelist[edge])
            setadjacent!(edgelist[Edge(edge[1],dots[j+4])],1+s,dots[5-j])
            setadjacent!(edgelist[Edge(dots[j+4],edge[2])],1+s,dots[3+j])
        end
        set!(edgelist, Edge(dots[SVector(3, 4)]), EdgeAttributes(I,degs[4], zero(P), false))
        set!(edgelist, Edge(dots[SVector(4, 5)]), EdgeAttributes(I,degs[5], zero(P), false))
        set!(trilist, Triangle(dots[SVector(1, 4, 5)]), TriangleAttributes{F,P}())
        set!(trilist, Triangle(dots[SVector(4, 2, 3)]), TriangleAttributes{F,P}())
        set!(trilist, Triangle(dots[SVector(4, 3, 5)]), TriangleAttributes{F,P}())
    end
    return nothing
end


"""
  $(SIGNATURES)

performs the refinement of Green marked triangles.   
"""
function refine_green!(t::Triangle{I}, mesh::HPMesh{F, I, P}, refaux::RefineAux{I, P}) where {F <: AbstractFloat, I <: Integer, P <: Integer}
    (; points, edgelist, trilist) = mesh
    (; i, degs, dots) = refaux
    dots[1:3] .= t
    edge = longestedge(t)
    degs[1:3] .= (degree(edgelist[e]) for e in edges(t))
    degs[4] = max(maximum(abs, degs[SVector(1, 3)] - degs[SVector(2, 1)]), one(P))
    ep = edgelist[edge]
    middle = seen(ep)
    s = middle>0
    if s
        dots[4] = middle
        seen!(ep,zero(I))
    else
        points[i[]] = SVector(sum(points[edge]) / 2)
        dots[4] = i[]
        seen!(ep, i[])
        m = tag(ep)
        set!(edgelist, Edge(dots[SVector(1, 4)]), EdgeAttributes(I,degs[1], zero(P), false))
        set!(edgelist, Edge(dots[SVector(4, 2)]), EdgeAttributes(I,degs[1], zero(P), false))
        i[] += 1
    end
    setadjacent!(edgelist[Edge(edge[1],dots[4])],1+s,dots[3])
    setadjacent!(edgelist[Edge(dots[4],edge[2])],1+s,dots[3])
    set!(edgelist, Edge(dots[SVector(3, 4)]), EdgeAttributes(I,degs[4], zero(P), false))
    set!(trilist, Triangle(dots[SVector(1, 4, 3)]), TriangleAttributes{F,P}())
    set!(trilist, Triangle(dots[SVector(4, 2, 3)]), TriangleAttributes{F,P}())
    nothing
end

"""
    $(SIGNATURES)

performs the refinement process of `mesh`. It is assumed that the triangles of `mesh` had already been marked. If there is no marked triangle, this functions does nothing.
"""
function refine!(mesh::HPMesh{F, I, P}) where {F <: AbstractFloat, I <: Integer, P <: Integer}
    (; points, edgelist, trilist) = mesh
    i = I(length(points) + 1)
    n_edgelist = count(ismarked, edgelist)
    append!(points, Vector{SVector{2, F}}(undef, n_edgelist))
    refaux = RefineAux(i, mesh)
    for t in triangles(trilist)
        if isred(trilist[t])
            refine_red!(t, mesh, refaux)
        elseif isblue(trilist[t])
            refine_blue!(t, mesh, refaux)
        elseif isgreen(trilist[t])
            refine_green!(t, mesh, refaux)
        end
    end
    filter!(!ismarked, mesh.trilist)
    filter!(!ismarked, mesh.edgelist)
    upgen!(mesh)
    nothing
end

function refine(mesh::HPMesh{F,I,P}) where {F,I,P}
    out = copy(mesh)
    refine!(out)
    out
end 

"""
    $(SIGNATURES)

checks if the values in `pt` satisfy the p_conformity condition: `p₁+p₂>=p₃`.
"""
function check_p_conformity(pt)
    return sum(pt) ≥ 2maximum(pt)
end

function check_p_conformity(m::HPMesh)
    for t in triangles(m)
        degs = degrees(t, m)
        if !check_p_conformity(degs)
            return false
        end
    end
    return true
end
"""
    $(SIGNATURES)

recursively checks the p_conformity of the triangles in `mesh`, incrementing the degrees when necessary. 
"""
function p_conformity!(mesh::HPMesh{F, I, P}) where {F, I, P}
    (; trilist) = mesh
    still = true
    while still
        still = false
        for t in triangles(trilist)
            if !p_conformity!(mesh, t, 10)
                still = true
            end
        end
    end
    return
end
#fallback case for boundary edges, where neighbor == nothing
p_conformity!(::HPMesh, ::Nothing, d) = true
function p_conformity!(mesh::HPMesh{F, I, P}, t::Triangle{I}, d) where {F, I, P}
    (; edgelist) = mesh
    p, eds = psortededges(t, mesh)
    out = false
    if check_p_conformity(p)
        out = true
    else
        if d > 0
            setdegree!(edgelist[eds[1]], p[3] - p[2])
            t₁ = tag(edgelist[eds[1]]) > 0 ? nothing : neighbor(mesh,t,eds[1])
            if p_conformity!(mesh, t₁, d - 1)
                out = true
            else
                setdegree!(edgelist[eds[1]], p[1])
                setdegree!(edgelist[eds[2]], p[3] - p[1])
                t₂ = tag(edgelist[eds[2]]) > 0 ? nothing : neighbor(mesh, t, eds[2])
                    if p_conformity!(mesh, t₂, d - 1)
                        out = true
                    else
                        setdegree!(edgelist[eds[2]], p[2])
                    end
            end
        end
    end
    return out
end

"""
    $(SIGNATURES)
    
Returns the neighbor of triangle `t` on the other side of edge `e`. Notice that this function should run ONLY on interior edges
"""
function neighbor(mesh::HPMesh{F, I, P}, t::Triangle{I}, e::Edge{I}) where {F, I, P}
    (;trilist,edgelist) = mesh
    ep = edgelist[e]
    tag(ep)>0 && error(ArgumentError("neighbor only works on interior edges. A boundary edge was passed."))
    adj = adjacents(ep)
    v = adj[1] ∈ t ? adj[2] : adj[1]
    tt = Triangle(e[1],e[2],v)
    _, u = gettoken(trilist, tt)
    gettokenvalue(keys(trilist), u)
end

"""
    setdegrees!(p::Function,mesh::HPMesh)
sets the degrees of the edges according to function `p`. `p` should evaluate on points on the domain and return values greater than 1.
"""
function setdegrees!(p::Function,mesh::HPMesh)
    (;points,edgelist) = mesh
    done = false
    while !done
        done = true
        for e in keys(edgelist)
            deg = degree(edgelist[e])
            med = (sum(points[e])/2)
            if deg < round(Int,p(med))
                setdegree!(edgelist[e],deg+1)
                done = false
            elseif deg > round(Int,p(med))
                setdegree!(edgelist[e],deg-1)
                done = false
            end
        end
        p_conformity!(mesh)
    end
end
###########################################################################################
### ESTIM Functions
"""
    $(SIGNATURES)

`estim` function for grading a mesh towards a set `S`, where `d` is the distance function to `S`. The resulting mesh grading is such that `hₜ ∼ h^(1/μ)` if `d(v)=0` for some `v` vertex of `t`, and `hₜ ∼ α*h*dₜ^(1-μ)` in other case, where `h`,`μ` and `α` are settable parameters and `d(t)` is the distance from `t` to `S`. `tol` is the tolerance for checking if `d(v)==0`.

This function should be passed to `mark!` for marking the triangles. 
"""
function estim_distance(vert; h = 0.2, μ = 1, α = 1.5, tol = 1.0e-12, dist::Function)
    ℓ = minimum(sum(abs2, vert[:, SVector(1, 2)] - vert[:, SVector(3, 3)], dims = 1)) |> sqrt
    d = maximum(d(v) for v in eachcol(vert))
    return ℓ > (d > tol ? h^(1 / μ) : α * h * d^(1 - μ))
end


"""
    $(SIGNATURES)

`estim` function for grading a mesh towards a point `pp`. The resulting mesh grading is such that for each triangle `t`: `hₜ ∼ h^(1/μ)` if `pp ∈ t`  and `hₜ ∼ α*h*d(t,pp)^(1-μ)` in other case. `h`,`μ` and `α` are settable parameters and `d(t,pp)` is the distance from `t` to `pp`.

This function should be passed to `mark!` for marking the triangles.

Note that `estim_origin(vert)` is a slightly more efficient version of this function when `pp` is the origin. 
"""
function estim_point(vert; h = 0.2, μ = 1, α = 1.5, center)
    ℓ = minimum(sum(abs2, vert[:, SVector(1, 2)] - vert[:, SVector(3, 3)], dims = 1)) |> sqrt
    d = maximum(sum(abs2, vert .- center, dims = 1)) |> sqrt
    if intriangle(center, vert) >= 0
        return ℓ > h^(1 / μ)
    else
        return ℓ > α * h * d^(1 - μ)
    end
end

"""
    $(SIGNATURES)

`estim` function for grading a mesh towards the origin. See the docs for `estim_point`.
"""
estim_origin(vert; args...) = estim_point(vert; center = SVector(0.0, 0.0), args...)


