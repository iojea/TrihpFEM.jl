"""Sol PLOT"""

@recipe PlotSol (sol::FESolution,) begin
    linewidth = 0.5
    annotate = false
    title = ""
end


function Makie.plot!(p::PlotSol)
    lift(p[1]) do sol
        (;mesh,vals) = sol
        (; points, trilist, edgelist) = mesh
        tris = hcat([Vector(t) for t in triangles(trilist)]...)'
        poly!(p, points,tris,color = vals,colormap=:coolwarm)
    end
    # hidedecorations!(p)
    return p
end
