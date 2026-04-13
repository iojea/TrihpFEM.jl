const AuxMeshData{P,F} = Dictionary{NTuple{3,P},AuxDegData{F}} where {P,F}

struct Measure{T<:HPTriangulation,A<:AuxMeshData,Q<:Quadrature}
    mesh::T
    aux::A
    sch::Q
end

function Measure(mesh::HPMesh{F,I,P},degsch) where {F,I,P}
    aux = AuxMeshData{P,F}()
    for t in triangles(mesh)
        d = degrees(t,mesh)
        isin,(_,_) = gettoken(aux,d)
        if !isin
            set!(aux,d,AuxDegData(F,d))
        end
    end
    sch = quadrature(F,P,Val(2),Val(degsch))
    Measure(mesh,aux,sch)
end

function Measure(mesh::HPMesh{F,I,P}) where {F,I,P}
    (;edgelist) = mesh
    maxdeg = maximum(Meshes.degree.(edgelist))
    Measure(mesh,2maxdeg+1)
end

mesh(m::Measure) = m.mesh
Meshes.domainmesh(m::Measure) = domainmesh(m.mesh)

function dof(ed::Edge,m::Measure)
    indofs,token = gettoken(m.mesh.dofs.by_edge,ed)
    indofs || throw(ArgumentError("No degrees of freedom are associated to the edge $ed.l Maybe you need to compute the degrees of freedom of the mesh, using `degrees_of_freedom!`."))
    return gettokenvalue(m.mesh.dofs.by_edge,token)
end
     
function dof(t::Triangle,m::Measure)
    indofs,token = gettoken(m.mesh.dofs.by_tri,t)
    indofs || throw(ArgumentError("No degrees of freedom are associated to the triangle $t.l Maybe you need to compute the degrees of freedom of the mesh, using `degrees_of_freedom!`."))
    return gettokenvalue(m.mesh.dofs.by_tri,token)
end

elements(m::Measure) = elements(m.mesh)
