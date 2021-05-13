using Revise
using FieldDstar

dstar = DstarSearch(Index(1, 1), Index(9, 9), PointGrid(10, 10))
open_list = dstar.open_list
compute_shortest_path(dstar)

