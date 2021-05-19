function neighbouring_index_pairs(pg::PointGrid, idx::Index)
    x_min = max(idx.x-1, 1)
    x_max = min(idx.x+1, pg.w)
    y_min = max(idx.y-1, 1)
    y_max = min(idx.y+1, pg.h)

    function inRange(idx::Index)
        idx.x < x_min && (return false)
        idx.x > x_max && (return false)
        idx.y < y_min && (return false)
        idx.y > y_max && (return false)
        return true
    end

    s1 = Index(idx.x + 1, idx.y + 0)
    s2 = Index(idx.x + 1, idx.y + 1)
    s3 = Index(idx.x + 0, idx.y + 1)
    s4 = Index(idx.x - 1, idx.y + 1)
    s5 = Index(idx.x - 1, idx.y + 0)
    s6 = Index(idx.x - 1, idx.y - 1)
    s7 = Index(idx.x + 0, idx.y - 1)
    s8 = Index(idx.x + 1, idx.y - 1)

    S = [s1, s2, s3, s4, s5, s6, s7, s8]
    valid_bools = [inRange(idx) & !pg.pred_coll(idx) for idx in S]
    is_valid_pair = (i, j) -> (valid_bools[i] & valid_bools[j])

    Channel() do c
        for i in 1:7
            is_valid_pair(i, i+1) && put!(c, (S[i], S[i+1]))
        end
        is_valid_pair(7, 1) && put!(c, (S[7], S[1]))
    end
end

function neighbouring_node_pairs(dstar::DstarSearch, node::Node)
    Channel() do c
        for idx_pair in neighbouring_index_pairs(dstar.pgrid, node.idx)
            s1 = get_node(dstar, idx_pair[1])
            s2 = get_node(dstar, idx_pair[2])
            put!(c, (s1, s2))
        end
    end
end
