using Revise
using FieldDstar
using GR

dstar = DstarSearch(Index(1, 1), Index(27, 27), PointGrid(30, 30))
open_list = dstar.open_list
compute_shortest_path(dstar, debug_visual=true)
