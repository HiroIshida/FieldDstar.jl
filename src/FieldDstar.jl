module FieldDstar

using StaticArrays
using DataStructures
using LinearAlgebra
using Plots
gr()

export PointGrid, Node, DstarSearch, Index, compute_shortest_path
export get_node, neighbouring_nodes, visualize, extract_path


macro debugassert(test)
  esc(:(if $(@__MODULE__).debugging()
    @assert($test)
   end))
end
debugging() = false

struct Index{T}
    x::T
    y::T
end
Index(x, y) = Index{Int}(x, y)
convert(::Type{Index{Float64}}, idx::Index{<:Int}) = Index{Float64}(idx.x, idx.y)

mutable struct Node
    g::Float64
    rhs::Float64
    idx::Index{Int}
    is_visited::Bool
end
Node(idx::Index{Int}) = Node(idx, false)
Node(idx::Index{Int}, is_visited::Bool) = Node(Inf, 0, idx, is_visited)

function euclidean_distance(node1::Node, node2::Node)
    x1, y1 = node1.idx.x, node1.idx.y
    x2, y2 = node2.idx.x, node2.idx.y
    return sqrt((x1 - x2)^2 + (y1 - y2)^2)
end

abstract type Heuristic end
struct EuclideanHeuristic <: Heuristic end
function (h::EuclideanHeuristic)(node1::Node, node2::Node) 
    return euclidean_distance(node1, node2)
end

struct KeyVal
    first::Float64
    second::Float64
end

function KeyVal(node_here::Node, node_start::Node, h::Heuristic)
    first = min(node_here.g, node_here.rhs) + h(node_here, node_start)
    second = min(node_here.g, node_here.rhs)
    return KeyVal(first, second)
end

function Base.isless(x::KeyVal, y::KeyVal)
    # see end of p. 10 of the paper
    x.first < y.first && (return true)
    (x.first == y.first && x.second < y.second) && (return true)
    return false
end

# unit cell grid
mutable struct PointGrid
    w::Int
    h::Int
    pred_coll::Function
end

function PointGrid(w::Int, h::Int, pred_coll::Function)
    return PointGrid(nodes, w, h, pred_coll)
end

mutable struct DstarSearch 
    s_start::Node
    s_goal::Node
    open_list::PriorityQueue{Node, KeyVal}
    nodes::Vector{Vector{Node}}
    pgrid::PointGrid
    is_field::Bool
end
@inline get_node(dstar::DstarSearch, idx::Index{Int}) = dstar.nodes[idx.x][idx.y]

function DstarSearch(idx_start::Index{Int}, idx_goal::Index{Int}, pgrid::PointGrid; is_field=true)
    nodes = [[Node(Index(i, j)) for j in 1:pgrid.h] for i in 1:pgrid.w]
    open_list = PriorityQueue{Node, KeyVal}()
    s_start = nodes[idx_start.x][idx_start.y]
    s_goal = nodes[idx_goal.x][idx_goal.y]
    return DstarSearch(s_start, s_goal, open_list, nodes, pgrid, is_field)
end

function DataStructures.enqueue!(dstar::DstarSearch, key::Node, val::KeyVal)
    key.is_visited = true;
    enqueue!(dstar.open_list, key, val)
end

function neighbouring_nodes(dstar::DstarSearch, node::Node)
    idx = node.idx
    pg = dstar.pgrid
    @debugassert idx.x > 0
    @debugassert idx.y > 0
    @debugassert idx.x <= pg.w
    @debugassert idx.y <= pg.h

    vec = Vector{Node}(undef, 0)
    for i in max(idx.x-1, 1):min(idx.x+1, pg.w)
        for j in max(idx.y-1, 1):min(idx.y+1, pg.h)
            (i==idx.x && j==idx.y) && continue
            pg.pred_coll(Index(i, j)) && continue
            node = get_node(dstar, Index{Int}(i, j))
            push!(vec, node)
        end
    end
    return vec
end

function compute_shortest_path(dstar::DstarSearch; debug_visual=false)
    function predicate_continue()
        s_start = dstar.s_start

        s_start.g != s_start.rhs && (return true)
        return false
    end

    itr = 0

    key = dstar.s_goal
    val = KeyVal(key, dstar.s_start, EuclideanHeuristic())
    enqueue!(dstar, key, val)

    while predicate_continue()
        s = dequeue!(dstar.open_list)
        if s.g > s.rhs
            s.g = s.rhs
            for s_ in neighbouring_nodes(dstar, s)
                update_state(dstar, s_)
            end
        else
            s.g = Inf
            (update_state(s_) for s_ in neighbouring_nodes(dstar, s))
            update_state(dstar, s)
        end
        debug_visual && visualize(dstar, s)
        itr += 1
    end
    return itr
end

function update_state(dstar::DstarSearch, s::Node)
    !s.is_visited && (s.g = Inf)
    if s!= dstar.s_goal
        nbrs = neighbouring_nodes(dstar, s)
        f(s_::Node) = euclidean_distance(s, s_) + s_.g
        if dstar.is_field
            s.rhs = minimum(compute_pair_ahead_rhs(pair)[1] for pair in neighbouring_node_pairs(dstar, s))
        else
            s.rhs = minimum(f(s_) for s_ in nbrs)
        end
    end
    haskey(dstar.open_list, s) && delete!(dstar.open_list, s)
    s.g != s.rhs && enqueue!(dstar, s, KeyVal(s, dstar.s_start, EuclideanHeuristic()))
end

function extract_path(dstar::DstarSearch)
    if dstar.is_field
        function backtrack(idx::Index{Int})
            node = get_node(dstar, idx)
            cost_min = Inf
            float_idx_best = nothing
            for pair in neighbouring_node_pairs(dstar, node)
                cost, float_idx = compute_pair_ahead_rhs(pair)
                if cost < cost_min
                    float_idx_best = float_idx
                    cost_min = cost
                end
            end
            return float_idx_best
        end

        function backtrack(idx::Index{Float64})
            error("not supported")
        end

        path = Vector{Index}(undef, 0)
        idx = dstar.s_start.idx
        push!(path, idx)
        while idx!=dstar.s_goal.idx
            idx = backtrack(idx)
        end
    else
        path = Vector{Index{Int}}(undef, 0)
        s = dstar.s_start
        push!(path, s.idx)
        while s!=dstar.s_goal
            nodes = collect(neighbouring_nodes(dstar, s))
            s = nodes[argmin([s.g for s in nodes])]
            push!(path, s.idx)
        end
        return path
    end
end

function visualize(dstar::DstarSearch, s::Node)
    gs = []
    for i in 1:dstar.pgrid.w
        for j in 1:dstar.pgrid.h
            if dstar.nodes[i][j].g == Inf
                push!(gs, 0.0)
            else
                push!(gs, dstar.nodes[i][j].g)
            end

            if i==s.idx.x && j==s.idx.y
                gs[end] = 0.0
            end
        end
    end
    xs = 1:dstar.pgrid.w
    ys = 1:dstar.pgrid.h
    GR.heatmap(xs, ys, reshape(gs, dstar.pgrid.w, dstar.pgrid.h))
end

include("piece_wise.jl")

end
