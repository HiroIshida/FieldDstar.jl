using Revise
using FieldDstar
using GR

dstar = DstarSearch(Index(1, 1), Index(9, 9), PointGrid(10, 10))
open_list = dstar.open_list
compute_shortest_path(dstar)

gs = []
for i in 1:10
    for j in 1:10
        if dstar.nodes[i][j].g == Inf
            push!(gs, 0.0)
        else
            push!(gs, dstar.nodes[i][j].g)
        end
    end
end

GR.heatmap(reshape(gs, 10, 10))


