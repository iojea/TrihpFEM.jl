"""
    $(SIGNATURES)

builds a `HPMesh` with elements of size `h` of the circle with center `center` and radius `rad`. If the center is ommited, it is assumed to be the origin. If the radius is also ommited, it is assumed to be 1. 
"""
function circmesh(center,rad,h)
    cx,cy = center
    n     = Int(1+2π*rad÷h)
    θ     = range(start=0.,stop=2π,length=n)[1:end-1]
    verts = [cx .+ cos.(θ');cy .+ sin.(θ')]
    hpmesh(verts,h)
end
circmesh(rad,h) = circmesh((0.,0.),rad,h)
circmesh(h) = circmesh(1.,h)

"""
   correct_boundary_circular(mesh::HPMesh)
corrects the boundary nodes of `mesh` projecting them to the boundary.  
"""
function correct_boundary_circular(mesh::HPMesh)
    (;points,edgelist) = mesh
    for e in keys(edgelist)
        if marker(edgelist[e])==1
            i,j = e
            # points[:,i] .= points[:,i]/norm(points[:,i])
            # points[:,j] .= points[:,j]/norm(points[:,j])
            points[i] = points[i]/norm(points[i])
            points[j] = points[j]/norm(points[j])
        end
    end
end

"""
    circmesh_graded_center(h,μ;maxiter,rec)
builds a mesh of a unitary disc, with triangles of edgelength `h`, and graded towards the center with graduation parameter `μ`.  
"""
function circmesh_graded_center(h,μ;maxiter=4,rec=false)
    k = 0
    mesh = circmesh(h)
    flag = true
    mshs = [copy(mesh)]
    while k≤maxiter && flag
        mark!(estim_origin,mesh,h=h,μ=μ)
        rec ? push!(mshs,copy(mesh)) : nothing
        if count(ismarked,mesh.trilist) > 0
            refine!(mesh)
            correct_boundary_circular(mesh)
            k += 1
            rec ? push!(mshs,copy(mesh)) : nothing
        else
            flag = false
        end
    end
    return rec ? mshs : mesh
end


"""
   l_mesh(h)
builds a mesh of an L shaped domain with triangles of edge-length `h`. 
"""
function l_mesh(h)
    x = range(start=0,stop=2,length=Int(1+2÷h))
    x0 = range(start=0,stop=1,length=Int(1+1÷h))    
    n = length(x)
    n0 = length(x0)
    y = range(start=0,stop=2,length=Int(1+2÷h))   
    y0 = range(start=0,stop=1,length=Int(1+1÷h))        
    m = length(y) 
    m0 = length(y0) 
    points   =  [-1. -1;1. -1.;1. 0.;0. 0.;0. 1;-1 1]'
    
    edgelist = hcat([[i,i+1] for i in 1:size(points,2)-1]...,
                          [size(points,2),1])
    marks    = Vector{Int8}(ones(Int8,size(edgelist,2)))
    tri = TriangulateIO(;pointlist=points,segmentlist=edgelist,segmentmarkerlist=marks)
    maxa = Printf.@sprintf "%.15f" h^2/2
    angle= Printf.@sprintf "%.15f" 30.
    (tri,_) = triangulate("ea$(maxa)q$(angle)pQ",tri)
    HPMesh(tri)
end


"""
   l_graded(h,μ;maxiter=6,rec=false)
builds a mesh of an L shaped domain with triangles of edge-length `h`, graded towards the inner vertex with graduation parameter `μ`.
"""
function l_graded(h,μ;maxiter=6,rec=false)
    k = 0
    mesh = l_mesh(h)
    flag = true
    mshs = [copy(mesh)]
    while k≤maxiter && flag
        mark!(estim_distance_origin_v,mesh,h=h,μ=μ)
        rec ? push!(mshs,copy(mesh)) : nothing
        if count(ismarked,mesh.trilist) > 0
            refine!(mesh)
            k += 1
            rec ? push!(mshs,copy(mesh)) : nothing
        else
            flag = false
        end
    end
    return rec ? mshs : mesh
end

