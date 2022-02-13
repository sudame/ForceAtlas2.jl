using Graphs, GraphPlot, ForceAtlas2

g = barabasi_albert(100, 3)
locs_x, locs_y = forceatlas2_layout(g, 1000)

gplot(g, locs_x, locs_y)