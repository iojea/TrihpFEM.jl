"""
$(SIGNATURES)

Is a cache for the refining process. It stores data that will be updated for each triangle.
"""
struct RefineAux{I<:Integer,P<:Integer}
    i::Base.RefValue{I}
    degs::MVector{6,P}
    dots::MVector{6,I}
    seen::Dictionary{HPEdge{I},I}
end

function RefineAux{I,P}() where {I<:Integer,P<:Integer}
    RefineAux(MVector{6,P}(zeros(6)),MVector{6,I}(zeros(6)),Dictionary{HPEdge{I},I}())
end

function RefineAux(i,mesh::HPMesh{F,I,P}) where {F<:AbstractFloat,I<:Integer,P<:Integer}
    seen = fill(zero(I),filter(ismarked,mesh.edgelist))
    RefineAux(Ref(i),MVector{6,P}(zeros(6)),MVector{6,I}(zeros(6)),seen)
end

"""
    $(SIGNATURES)

Marks triangles with vertices `vert` for which `estim(vert,estim_param)` returns `true`, and then run the `_h_conformity!` routine in order to propagate the markings to adjacent triangles as needed. 
"""
function mark!(estim,mesh::HPMesh{F,I,P};estim_params...) where {F,I,P}
    (;points,trilist,edgelist) = mesh
    tri  = MMatrix{2,3}(zero(F) for i in 1:2,j in 1:3)
    for t in triangles(trilist)
        for (i,tt) in enumerate(t)
            tri[:,i] .= points[tt]
        end
        if estim(tri;estim_params...)
            mark!.(getindices(edgelist,edges(t))) #AQuí había un FOR que cambié
        end  
    end
    _h_conformity!(mesh)
end

"""
    $(SIGNATURES)

checks if a point `p` belongs to thte triangle with vertices `a`,`b` and `c`. 
"""
function intriangle(p::T,a::V,b::V,c::V) where {T<:AbstractArray,V<:AbstractArray}
    if abs(orient(a,b,p)+orient(b,c,p)+orient(c,a,p))==3
        return 1
    elseif orient(a,b,p)*orient(b,c,p)*orient(c,a,p) == 0
        return 0
    else
        return -1
    end
end
intriangle(p::T,vert::M) where {T<:AbstractArray,M<:AbstractArray} = intriangle(p,eachcol(vert)...)



"""
    $(SIGNATURES)

Performs the marking of previously un-marked triangles in order to avoid hanging nodes.
"""

function _h_conformity!(mesh::HPMesh)
    (;trilist,edgelist) = mesh
    still = true
    while still
        still = false
        for t in triangles(trilist)
            num_marked = count(ismarked,getindices(edgelist,edges(t)))     
            if num_marked>0
                long_edge = edgelist[longestedge(t)]
                if !ismarked(long_edge)
                    mark!(long_edge)
                    still = true
                    num_marked += 1
                end
                mark!(trilist[t],num_marked)
            end
        end
    end
end


"""
  $(SIGNATURES)

performs the refinement of Red marked triangles.   
"""
function refine_red!(t::HPTriangle{I},mesh::HPMesh{F,I,P},refaux::RefineAux{I,P}) where {F<:AbstractFloat,I<:Integer,P<:Integer}
    (;points,edgelist,trilist) = mesh
    (;i,degs,dots,seen) = refaux
    dots[1:3] .= t
    t_edges    = edges(t) 
    degs[1:3] .= (degree(edgelist[e]) for e in t_edges)
    degs[4:6] .= max.(abs.(degs[SVector(1,2,3)]-degs[SVector(3,1,2)]),one(P))
    for j in eachindex(t_edges)
        edge = t_edges[j]
        k    = seen[edge]
        if k>0
            dots[j+3] = k
        else
            points[i[]] = HPPoint(sum(points[edge])/2)
            dots[j+3]   = i[]
            set!(seen,edge,i[])
            m = marker(edgelist[edge])
            set!(edgelist,HPEdge(edge[1],i[]),EdgeProperties(degs[j],m,false))
            set!(edgelist,HPEdge(i[],edge[2]),EdgeProperties(degs[j],m,false))
            i[] += 1
        end
    end 
    set!(edgelist,HPEdge(dots[SVector(6,4)]),EdgeProperties(degs[4],zero(P),false))
    set!(edgelist,HPEdge(dots[SVector(4,5)]),EdgeProperties(degs[5],zero(P),false))    
    set!(edgelist,HPEdge(dots[SVector(5,6)]),EdgeProperties(degs[6],zero(P),false))
    set!(trilist,HPTriangle(dots[SVector(1,4,6)]),TriangleProperties(P,F))
    set!(trilist,HPTriangle(dots[SVector(4,2,5)]),TriangleProperties(P,F))
    set!(trilist,HPTriangle(dots[SVector(6,5,3)]),TriangleProperties(P,F))
    set!(trilist,HPTriangle(dots[SVector(5,6,4)]),TriangleProperties(P,F))
end

"""
  $(SIGNATURES)

performs the refinement of Blue marked triangles.   
"""
function refine_blue!(t::HPTriangle{I},mesh::HPMesh{F,I,P},refaux::RefineAux{I,P}) where {F<:AbstractFloat,I<:Integer,P<:Integer}
    (;points,edgelist,trilist) = mesh
    (;i,degs,dots,seen) = refaux
    dots[1:3] .= t
    t_edges    = edges(t)
    degs[1:3] .= (degree(edgelist[e]) for e in t_edges)
    degs[4]    = max(maximum(abs,degs[SVector(1,3)]-degs[SVector(2,1)]),one(P))
    if ismarked(edgelist[t_edges[2]])
        degs[5] = max(maximum(abs,degs[SVector(1,2)]-degs[SVector(2,4)]),one(P))
        for j in 1:2
            edge = t_edges[j]
            k    = seen[edge]
            if k>0
                dots[j+3] = k
            else
                points[i[]] = HPPoint(sum(points[edge])/2)
                dots[j+3]   = i[]
                set!(seen,edge,i[])
                m = marker(edgelist[edge])
                set!(edgelist,HPEdge(edge[1],i[]),EdgeProperties(degs[j],m,false))
                set!(edgelist,HPEdge(i[],edge[2]),EdgeProperties(degs[j],m,false))
                i[] += 1
            end
        end
        set!(edgelist,HPEdge(dots[SVector(3,4)]),EdgeProperties(degs[4],zero(P),false))
        set!(edgelist,HPEdge(dots[SVector(5,4)]),EdgeProperties(degs[5],zero(P),false)) 
        set!(trilist,HPTriangle(dots[SVector(1,4,3)]),TriangleProperties(P,F))
        set!(trilist,HPTriangle(dots[SVector(4,2,5)]),TriangleProperties(P,F))
        set!(trilist,HPTriangle(dots[SVector(4,5,3)]),TriangleProperties(P,F))
    elseif ismarked(edgelist[t_edges[3]])
        degs[5]    = max(maximum(abs,degs[SVector(1,3)]-degs[SVector(3,4)]),one(P))
        for j in 0:1
            edge = t_edges[1+2j]
            k    = seen[edge]
            if k>0
                dots[j+4] = k
            else
                # points[i[]] = HPPoint(sum(points[edge])/2.)
                points[i[]] = HPPoint(sum(points[edge])/2)
                dots[j+4]   = i[]
                set!(seen,edge,i[])
                m = marker(edgelist[edge])
                set!(edgelist,HPEdge(edge[1],i[]),EdgeProperties(degs[1+2j],m,false))
                set!(edgelist,HPEdge(i[],edge[2]),EdgeProperties(degs[1+2j],m,false))
                i[] += 1
            end
        end
        set!(edgelist,HPEdge(dots[SVector(3,4)]),EdgeProperties(degs[4],zero(P),false))
        set!(edgelist,HPEdge(dots[SVector(4,5)]),EdgeProperties(degs[5],zero(P),false))    
        set!(trilist,HPTriangle(dots[SVector(1,4,5)]),TriangleProperties(P,F))
        set!(trilist,HPTriangle(dots[SVector(4,2,3)]),TriangleProperties(P,F))
        set!(trilist,HPTriangle(dots[SVector(4,3,5)]),TriangleProperties(P,F))
    end
    return  nothing
end


"""
  $(SIGNATURES)

performs the refinement of Green marked triangles.   
"""
function refine_green!(t::HPTriangle{I},mesh::HPMesh{F,I,P},refaux::RefineAux{I,P}) where {F<:AbstractFloat,I<:Integer,P<:Integer}
    (;points,edgelist,trilist) = mesh
    (;i,degs,dots,seen) = refaux
    dots[1:3] .= t
    edge = longestedge(t)
    degs[1:3] .= (degree(edgelist[e]) for e in edges(t))
    degs[4]    = max(maximum(abs,degs[SVector(1,3)]-degs[SVector(2,1)]),one(P))
    k = seen[edge]
    if k>0
        dots[4] = k
    else
        points[i[]] = HPPoint(sum(points[edge])/2)
        dots[4]     = i[]
        set!(seen,edge,i[])
        oldedge = edgelist[edge] 
        m = marker(oldedge)
        set!(edgelist,HPEdge(dots[SVector(1,4)]),EdgeProperties(degs[1],zero(P),false))
        set!(edgelist,HPEdge(dots[SVector(4,2)]),EdgeProperties(degs[1],zero(P),false))
        i[] += 1
    end
    set!(edgelist,HPEdge(dots[SVector(3,4)]),EdgeProperties(degs[4],zero(P),false))
    set!(trilist,HPTriangle(dots[SVector(1,4,3)]),TriangleProperties(P,F))
    set!(trilist,HPTriangle(dots[SVector(4,2,3)]),TriangleProperties(P,F))
end

"""
    $(SIGNATURES)

performs the refinement process of `mesh`. It is assumed that the triangles of `mesh` had already been marked. If there is no marked triangle, this functions does nothing.
"""
function refine!(mesh::HPMesh{F,I,P}) where {F<:AbstractFloat,I<:Integer,P<:Integer}
    (;points,edgelist,trilist) = mesh
    i = I(length(points)+1)
    n_edgelist = count(ismarked,edgelist)
    #n_tris  = 4count(isred,trilist)+3count(isblue,trilist)+2count(isgreen,trilist)
    # sizehint!(trilist,length(trilist)+n_tris)
    # sizehint!(edgelist,length(trilist)+n_edgelist)
    append!(points,Vector{HPPoint{2,F}}(undef,n_edgelist))
    refaux  = RefineAux(i,mesh)
    for t in triangles(trilist)
        if isred(trilist[t])
            refine_red!(t,mesh,refaux)
        elseif isblue(trilist[t])
            refine_blue!(t,mesh,refaux)
        elseif isgreen(trilist[t])
            refine_green!(t,mesh,refaux)
        end
    end
    filter!(!ismarked,mesh.trilist)
    filter!(!ismarked,mesh.edgelist)
end


"""
    $(SIGNATURES)

checks if the values in `pt` satisfy the p_conformity condition: `p₁+p₂>=p₃`.
"""
function check_p_conformity(pt::Vector{T}) where T<:Integer
    sum(pt)  ≥ 2maximum(pt) 
end

"""
    $(SIGNATURES)

recursively checks the p_conformity of the triangles in `mesh`, incrementing the degrees when necessary. 
"""
function p_conformity!(mesh::HPMesh{F,I,P}) where {F,I,P}
    (;trilist) = mesh
    for t in triangles(trilist)
        p_conformity!(mesh,t,10)
    end
end
function p_conformity!(mesh::HPMesh{F,I,P},t::HPTriangle{I},d) where {F,I,P}
    (;edgelist) = mesh
    p,eds = pedges(t)
    out = false
    if check_p_conformity(p)
        out = true
    else
        if d>0
            setdegree!(edgelist[eds[1]],p[1] + 1)
            t₁ = neighbor(mesh,t,eds[1])
            if p_conformity!(mesh,t₁,d-1)
                out = true
            else
                setdegree!(edgelist[eds[1]],p[1])
                setdegree!(edgelist[eds[2]],p[2]+1)
                t₂ = neighbor(mesh,t,eds[2])
                if p_conformity!(mesh,t₂,d-1)
                    out = true
                else
                    setdegree!(edgelist[eds[2]],p[2])
                end
            end
        end
    end
    out
end

"""
    $(SIGNATURES)
    
If it exists, returns the neighbor of triangle `t` along the edge `e`. If `e` is a boundary edge, it returns `nothing`:w
.
"""
function neighbor(mesh::HPMesh{F,I,P},t::HPTriangle{I},e::HPEdge{I}) where {F,I,P}
    for tb in triangles(mesh)
        if  (e in edges(tb)) && t!=tb
            return tb
        end
    end
    return nothing
end

###########################################################################################
### ESTIM Functions
"""
    $(SIGNATURES)

`estim` function for grading a mesh towards a set `S`, where `d` is the distance function to `S`. The resulting mesh grading is such that `hₜ ∼ h^(1/μ)` if `d(v)=0` for some `v` vertex of `t`, and `hₜ ∼ α*h*dₜ^(1-μ)` in other case, where `h`,`μ` and `α` are settable parameters and `d(t)` is the distance from `t` to `S`. `tol` is the tolerance for checking if `d(v)==0`.

This function should be passed to `mark!` for marking the triangles. 
"""
function estim_distance(vert;h=0.2,μ=1,α=1.5,tol=1e-12,dist::Function)
    ℓ = minimum(sum(abs2,vert[:,SVector(1,2)]-vert[:,SVector(3,3)],dims=1)) |> sqrt
    d = maximum(d(v) for v in eachcol(vert))
    ℓ > (d>tol ? h^(1/μ) : α*h*d^(1-μ)) 
end


"""
    $(SIGNATURES)

`estim` function for grading a mesh towards a point `pp`. The resulting mesh grading is such that for each triangle `t`: `hₜ ∼ h^(1/μ)` if `pp ∈ t`  and `hₜ ∼ α*h*d(t,pp)^(1-μ)` in other case. `h`,`μ` and `α` are settable parameters and `d(t,pp)` is the distance from `t` to `pp`.

This function should be passed to `mark!` for marking the triangles.

Note that `estim_origin(vert)` is a slightly more efficient version of this function when `pp` is the origin. 
"""
function estim_point(vert;h=0.2,μ=1,α=1.5,center)
    ℓ = minimum(sum(abs2,vert[:,SVector(1,2)]-vert[:,SVector(3,3)],dims=1)) |> sqrt
    d = maximum(sum(abs2,vert.-center,dims=1)) |> sqrt
    if intriangle(center,vert) >= 0
        return ℓ > h^(1/μ)
    else
        return ℓ > α*h*d^(1-μ)
    end
end

"""
    $(SIGNATURES)

`estim` function for grading a mesh towards the origin. See the docs for `estim_point`.
"""
estim_origin(vert;args...) = estim_point(vert;center=ORIGIN,args...)


