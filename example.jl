using Revise
using FieldDstar
using GR

function pred_coll(idx::Index)
    if 10 < idx.x && idx.x < 20
        if 10 < idx.y && idx.y < 20
            return true
        end
    end
    return false
end

dstar = DstarSearch(Index(1, 1), Index(27, 27), PointGrid(30, 30, pred_coll))
open_list = dstar.open_list
compute_shortest_path(dstar, debug_visual=false)
path = extract_path(dstar)
all(!pred_coll(s.idx) for s in path)

