module FieldDstar

using StaticArrays
using DataStructures
using LinearAlgebra

export PointGrid, Node, DstarSearch, Index, compute_shortest_path


macro debugassert(test)
  esc(:(if $(@__MODULE__).debugging()
    @assert($test)
   end))
end
debugging() = false

struct Index
    x::Int
    y::Int
end

mutable struct Node
    g::Float64
    rhs::Float64
    idx::Index
end
Node(idx::Index) = Node(Inf, 0, idx)

function euclidean_distance(s1::Node, s2::Node)
    x1, y1 = s1.idx.x, s1.idx.y
    x2, y2 = s2.idx.x, s2.idx.y
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

struct KeyValOrdering <: Base.Order.Ordering end
import Base.Order.lt
function lt(o::KeyValOrdering, a::KeyVal, b::KeyVal)
    return true
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
end

function PointGrid(w::Int, h::Int)
    return PointGrid(nodes, w, h)
end

function neighbouring_indexes(pg::PointGrid, idx::Index)
    @debugassert idx[1] > 0
    @debugassert idx[2] > 0
    @debugassert idx[1] <= pg.w
    @debugassert idx[2] <= pg.h

    Channel() do c # like python's generator
        for i in max(idx[1], 1):min(idx[1], pg.w)
            for j in max(idx[2], 1):min(idx[2], pg.w)
                put!(c, Index(i, j))
            end
        end
    end
end

mutable struct DstarSearch 
    s_start::Node
    s_goal::Node
    open_list::PriorityQueue{Node, KeyVal}
    nodes::Vector{Vector{Node}}
    pgrid::PointGrid
end
@inline get_node(dstar::DstarSearch, idx::Index) = dstar.nodes[idx.x][idx.y]

function DstarSearch(idx_start::Index, idx_goal::Index, pgrid::PointGrid)
    nodes = [[Node(Index(i, j)) for j in 1:pgrid.h] for i in 1:pgrid.w]
    open_list = PriorityQueue{Node, KeyVal}(KeyValOrdering())
    s_start = nodes[idx_start.x][idx_start.y]
    s_goal = nodes[idx_goal.x][idx_goal.y]
    key = nodes[idx_goal.x][idx_goal.y]
    val = KeyVal(key, s_start, EuclideanHeuristic())
    enqueue!(open_list, key, val)
    return DstarSearch(s_start, s_goal, open_list, nodes, pgrid)
end

function neighbouring_nodes(dstar::DstarSearch, node::Node)
    Channel() do c # like a python generator
        for idx_neigh in neighbouring_indexes(dstar.pgrid, node.idx)
            put!(c, get_node(dstar, idx_neigh)A)
        end
    end
end

function compute_shortest_path(dstar::DstarSearch)
    function predicate_continue()
        s_start = dstar.s_start
        minval = minimum(pair.second for pair in dstar.open_list.xs)
        cond1 = minval < KeyVal(s_start, s_start, EuclideanHeuristic())
        cond2 = s_start.g == s_start.rhs
        return (cond1 || cond2)
    end

    while predicate_continue()
        s = dequeue!(dstar.open_list)
        if s.g > s.rhs
            s.g = s.rhs
            (update_state(n) for n in neighbouring_nodes())
        else
            s.g = Inf
            (update_state(n) for n in neighbouring_nodes())
            update_state(s)
        end
    end
end

function update_state(dstar::DstarSearch, s::Node)
    if s!= dstar.s_goal
        nbrs = neighbouring_nodes(dstar, s)
        f(s_::Node) = euclidean_distance(s, s_) + s_.g
        s.rhs = minimum(f(s_) for s_ in nbrs)
    end
    haskey(dstar.open_list, s) && delete(dstar.open_list, s)
    s.g != s.rhs && enqueue!(dstar.open_list, s)
end

end # module
