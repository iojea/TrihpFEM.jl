const AuxMeshData{P,F} = Dictionary{NTuple{3,P},AuxDegData{F}} where {P,F}

struct Measure{T<:HPTriangulation,A<:AuxMeshData,Q<:Quadrature}
    mesh::T
    aux::A
    sch::Q
end
Measure(mesh,aux,sch) = Mesh{typeof(mesh),typeof(aux),typeof(sch)}(mesh,aux,sch)

function Measure(mesh::HPMesh{F,I,P},degsch) where {F,I,P}
    aux = AuxMeshData{P,F}()
    for t in triangles(mesh)
        d = degrees(t,mesh)
        isin,(_,_) = gettoken(aux,d)
        if !isin
            set!(aux,d,AuxDegData(F,d))
        end
    end
    degsch = isodd(degsch) ? P(degsch) : P(degsch+1)
    sch = gmquadrature(F,Val(2),degsch)
    Measure(mesh,aux,sch)
end

function Measure(mesh::HPMesh{F,I,P}) where {F,I,P}
    (;edgelist) = mesh
    maxdeg = maximum(Meshes.degree.(edgelist))
    Measure(mesh,2maxdeg+1)
end

Meshes.domainmesh(m::Measure) = domainmesh(m.mesh)



