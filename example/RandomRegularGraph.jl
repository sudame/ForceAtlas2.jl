using Graphs, GraphPlot, ForceAtlas2
using Cairo # optional 

g = random_regular_graph(100, 2)
locs_x, locs_y = forceatlas2_layout(g, 1000)

gplot(g, locs_x, locs_y)
