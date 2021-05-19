using Test 
using FieldDstar

function pred_coll(idx::Index)
    if 10 < idx.x && idx.x < 20
        if 10 < idx.y && idx.y < 20
            return true
        end
    end
    return false
end
map = PointGrid(30, 30, pred_coll)
dstar = DstarSearch(Index(1, 1), Index(27, 27), map; is_field=false)
field_dstar = DstarSearch(Index(1, 1), Index(27, 27), map; is_field=true)

@testset "neighbouring pairs" begin
    pairs = collect(FieldDstar.neighbouring_node_pairs(dstar, get_node(dstar, Index(1, 1))))
    @test length(pairs) == 2

    pairs = collect(FieldDstar.neighbouring_node_pairs(dstar, get_node(dstar, Index(10, 10))))
    @test length(pairs) == 6

    pairs = collect(FieldDstar.neighbouring_node_pairs(dstar, get_node(dstar, Index(10, 12))))
    @test length(pairs) == 4

    pairs = collect(FieldDstar.neighbouring_node_pairs(dstar, get_node(dstar, Index(14, 14))))
    @test length(pairs) == 0
end

@testset "planning" begin
    for planner in  [dstar, field_dstar]
        open_list = planner.open_list
        compute_shortest_path(planner, debug_visual=false)
        path = extract_path(planner)
        @test path[1] == planner.s_start
        @test path[end] == planner.s_goal
        @test all(!pred_coll(s.idx) for s in path)
    end
end
