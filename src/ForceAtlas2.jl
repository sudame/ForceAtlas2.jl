module ForceAtlas2
using Graphs

include("./Utils.jl")
export Node, linear_repulsion!


function init(g::AbstractGraph)
  _nodes = Node[]
  for v in vertices(g)
    node = Node(
      mass = degree(g, v),
      old_dx = 0,
      old_dy = 0,
      dx = 0,
      dy = 0,
      x = rand(),
      y = rand(),
    )
    push!(_nodes, node)
  end

  _edges = Edge[]
  for e in edges(g)
    edge = Edge(
      node1 = e.src,
      node2 = e.dst,
      weight = 1.0, # TODO: support weighted graphs
    )
    push!(_edges, edge)
  end

  return (_nodes, _edges)
end


function forceatlas2(
  g::AbstractGraph,
  iterations::Int = 100;

  # behavior alternatives
  outbound_attraction_distribution::Bool = false,
  linear_log_mode::Bool = false, # TODO: support linear log mode
  adjust_sizes::Bool = false, # TODO: support adjust sizes
  edge_weight_influence::Real = 1.0, # MEMO: This property doesn't have any effect because we don't support wighted graph now. 

  # performance
  jitter_tolerance::Real = 1.0,
  barnes_hut_optimize::Bool = true,
  barnes_hut_theta::Real = 1.2,
  multi_threaded::Bool = false, # TODO: support multi threading

  # tuning
  scaling_ratio::Real = 2.0,
  strong_gravity_mode::Bool = false,
  gravity::Real = 1.0
)
  speed = 1.0
  speed_efficiency = 1.0
  nodes, edges = init(g)
  outbound_attraction_compensation = 1.0
  if outbound_attraction_distribution
    outbound_attraction_compensation = mean([n.mass for n in nodes])
  end

  for _ = 1:iterations
    for n in nodes
      n.old_dx = n.dx
      n.old_dy = n.dy
      n.dx = 0
      n.dy = 0
    end

    if barnes_hut_optimize
      rootRegion = Region(nodes)
      build_sub_regions!(rootRegion)
      apply_force_on_nodes!(rootRegion, nodes, barnes_hut_theta; coefficient = scaling_ratio)
    else
      apply_repulsion!(nodes, scaling_ratio)
    end

    apply_gravity!(nodes, gravity, scaling_ratio; use_strong_gravity = strong_gravity_mode)

    apply_attraction!(nodes, edges, outbound_attraction_distribution, outbound_attraction_compensation, edge_weight_influence)

    asaf_result::AdjustSpeedAndApplyForcesResult = adjust_speed_and_apply_forces!(nodes, speed, speed_efficiency, jitter_tolerance)
    speed = asaf_result.speed
    speed_efficiency = asaf_result.speed_efficiency
  end

  return ([n.x for n in nodes], [n.y for n in nodes])
end

end
