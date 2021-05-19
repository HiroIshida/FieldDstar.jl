struct NodePair 
    s_mid::Node
    s_tip::Node
end

function compute_pair_ahead_rhs(pair::NodePair)
    s1, s2 = pair.s_mid, pair.s_tip
    g1, g2 = s1.g, s2.g
    cost1 = g1 + 1.0
    cost2 = g2 + sqrt(2.0)

    c = 1.0

    f = g1 - g2
    if f<=0.0
        float_idx = convert(Index{Float64}, s1.idx)
        return c + g1, float_idx
    end

    if f < c
        t = min(f/sqrt(c^2 - f^2), 1.0) # y in the paper
        float_idx_x = s1.idx.x * t + s2.idx.x * (1. - t)
        float_idx_y = s1.idx.y * t + s2.idx.y * (1. - t)
        float_idx = Index{Float64}(float_idx_x, float_idx_y) 
        return c * sqrt(1 + t^2) + f * (1. - t) + g2, float_idx
    end

    float_idx = convert(Index{Float64}, s2.idx)
    return c * sqrt(2.0) + g2, float_idx
end

function neighbouring_node_pairs(dstar::DstarSearch, node::Node)
    pg = dstar.pgrid
    idx = node.idx

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

    function add_node_pair(vec::Vector{NodePair}, i, j)
        is_valid_pair = valid_bools[i] & valid_bools[j]
        if is_valid_pair
            node1 = get_node(dstar, S[i])
            node2 = get_node(dstar, S[j])
            push!(vec, NodePair(node1, node2))
        end
    end

    vec = Vector{NodePair}(undef, 0)

    for i in [3, 5, 7]
        add_node_pair(vec, i, i-1)
        add_node_pair(vec, i, i+1)
    end
    # corner case
    add_node_pair(vec, 1, 7)
    add_node_pair(vec, 1, 2)
    return vec
end

function neighbouring_node_pairs(dstar::DstarSearch, idx::Index{Float64})
    # TODO super inefficient!
    # duplicate 
    idxes = Index{Int}[]
    if idx.x - round(idx.x) == 0
        i = round(idx.x)
        j_low = floor(idx.y) 
        j_high = ceil(idx.y)
        idx1 = Index(i+0, j_low)
        idx2 = Index(i+1, j_low)
        idx3 = Index(i+1, j_high)
        idx4 = Index(i+0, j_high)
        idx5 = Index(i-1, j_high)
        idx6 = Index(i-1, j_low)
        idxes = [idx1, idx2, idx3, idx4, idx5, idx6]
    else
        j = round(idx.y)
        i_low = floor(idx.x) 
        i_high = ceil(idx.x)
        idx1 = Index(i_high, j+0)
        idx2 = Index(i_high, j+1)
        idx3 = Index(i_low, j+1)
        idx4 = Index(i_low, j+0)
        idx5 = Index(i_low, j-1)
        idx6 = Index(i_high, j-1)
        idxes = [idx1, idx2, idx3, idx4, idx5, idx6]
    end

    valid_bools = [inRange(idx) & !pg.pred_coll(idx) for idx in idxes]

    function add_node_pair(vec::Vector{NodePair}, i, j)
        is_valid_pair = valid_bools[i] & valid_bools[j]
        if is_valid_pair
            node1 = get_node(dstar, idxes[i])
            node2 = get_node(dstar, idxes[j])
            push!(vec, NodePair(node1, node2))
        end
    end

    vec = Vector{NodePair}(undef, 0)

    for i in [3, 5, 7]
        add_node_pair(vec, i, i-1)
        add_node_pair(vec, i, i+1)
    end
    # corner case
    add_node_pair(vec, 1, 7)
    add_node_pair(vec, 1, 2)
    return vec
end
