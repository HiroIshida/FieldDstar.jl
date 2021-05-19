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


@testset "neighbouring pairs" begin
    pairs = collect(FieldDstar.neighbouring_index_pairs(map, Index(1, 1)))
    @test length(pairs) == 2

    pairs = collect(FieldDstar.neighbouring_index_pairs(map, Index(10, 10)))
    @test length(pairs) == 6

    pairs = collect(FieldDstar.neighbouring_index_pairs(map, Index(10, 12)))
    @test length(pairs) == 4

    pairs = collect(FieldDstar.neighbouring_index_pairs(map, Index(14, 14)))
    @test length(pairs) == 0
end

@testset "planning" begin
    dstar = DstarSearch(Index(1, 1), Index(27, 27), map)
    open_list = dstar.open_list
    compute_shortest_path(dstar, debug_visual=false)
    path = extract_path(dstar)
    @test path[1] == dstar.s_start
    @test path[end] == dstar.s_goal
    @test all(!pred_coll(s.idx) for s in path)
end
