using Base

Base.@kwdef mutable struct Node
  mass::Real
  old_dx::Real
  old_dy::Real
  dx::Real
  dy::Real
  x::Real
  y::Real
end

Base.@kwdef mutable struct Region
  mass::Real
  mass_center_x::Real
  mass_center_y::Real
  size::Real
  nodes::Vector{Node}
  subregions::Vector{Region}

  function Region(nodes::Vector{Node})
    mass = sum(map(node -> node.mass, nodes))
    mass_sum_x = sum(map(node -> node.mass * node.x, nodes))
    mass_sum_y = sum(map(node -> node.mass * node.y, nodes))

    mass_center_x = mass_sum_x / mass
    mass_center_y = mass_sum_y / mass

    size = maximum(map(n -> sqrt((n.x - mass_center_x)^2 + (n.y - mass_center_y)^2), nodes))

    new(
      mass,
      mass_center_x,
      mass_center_y,
      size,
      nodes,
      Region[],
    )
  end
end

Base.@kwdef mutable struct Edge
  node1::Int
  node2::Int
  weight::Real
end

function linear_repulsion!(n1::Node, n2::Node; coefficient::Real = 0)
  xdist = n1.x - n2.x
  ydist = n1.y - n2.y
  distance2 = xdist^2 + ydist^2

  if distance2 > 0
    factor = coefficient * n1.mass * n2.mass / distance2
    n1.dx += xdist * factor
    n1.dy += ydist * factor
    n2.dx -= xdist * factor
    n2.dy -= ydist * factor
  end
end

function linear_repulsion_region!(n::Node, r::Region; coefficient::Real = 0)
  xdist = n.x - r.mass_center_x
  ydist = n.y - r.mass_center_y
  distance2 = xdist^2 + ydist^2

  if distance2 > 0
    factor = coefficient * n.mass * r.mass / distance2
    n.dx += xdist * factor
    n.dy += ydist * factor
  end
end

function linear_gravity!(n::Node, g::Real)
  xdist = n.x
  ydist = n.y
  distance = sqrt(xdist^2 + ydist^2)

  if distance > 0
    factor = n.mass * g / distance
    n.dx -= xdist * factor
    n.dy -= ydist * factor
  end
end

function strong_gravity!(n::Node, g::Real; coefficient::Real = 0)
  xdist = n.x
  ydist = n.y

  if xdist != 0 && ydist != 0
    factor = coefficient * n.mass * g
    n.dx -= xdist * factor
    n.dy -= ydist * factor
  end
end

function linear_attraction!(n1::Node, n2::Node, e::Real, distributed_attraction::Bool; coefficient::Real = 0)
  xdist = n1.x - n2.x
  ydist = n1.y - n2.y

  factor::Real = 0
  if !distributed_attraction
    factor = -coefficient * e
  else
    factor = -coefficient * e / n1.mass
  end

  n1.dx += xdist * factor
  n1.dy += ydist * factor
  n2.dx -= xdist * factor
  n2.dy -= ydist * factor
end

function apply_repulsion!(nodes::Vector{Node}, coefficient::Real)
  i = 0
  for n1 in nodes
    j = i
    for n2 in nodes
      if j == 0
        break
      end
      linear_repulsion!(n1, n2, coefficient = coefficient)
      j -= 1
    end
    i += 1
  end
end

function apply_gravity!(nodes::Vector{Node}, gravity::Real, scaling_ratio::Real; use_strong_gravity::Bool = false)
  if !use_strong_gravity
    for n in nodes
      linear_gravity!(n, gravity)
    end
  else
    for n in nodes
      strong_gravity!(n, gravity, coefficient = scaling_ratio)
    end
  end
end

function apply_attraction!(nodes::Vector{Node}, edges::Vector{Edge}, distributed_attraction::Bool, coefficient::Real, edge_weight_influence::Real)
  if edge_weight_influence == 0
    for edge in edges
      linear_attraction!(
        nodes[edge.node1],
        nodes[edge.node2],
        1,
        distributed_attraction;
        coefficient = coefficient
      )
    end
  elseif edge_weight_influence == 1
    for edge in edges
      linear_attraction!(nodes[edge.node1],
        nodes[edge.node2],
        edge.weight,
        distributed_attraction;
        coefficient = coefficient
      )
    end
  else
    for edge in edges
      linear_attraction!(nodes[edge.node1],
        nodes[edge.node2],
        edge.weight^edge_weight_influence,
        distributed_attraction;
        coefficient = coefficient
      )
    end
  end
end

function build_sub_regions!(r::Region)
  if length(r.nodes) > 1
    topleft_nodes::Vector{Node} = Node[]
    bottomleft_nodes::Vector{Node} = Node[]
    topright_nodes::Vector{Node} = Node[]
    bottomright_nodes::Vector{Node} = Node[]

    for n in r.nodes
      if n.x < r.mass_center_x
        if n.y < r.mass_center_y
          push!(bottomleft_nodes, n)
        else
          push!(topleft_nodes, n)
        end
      else
        if n.y < r.mass_center_y
          push!(bottomright_nodes, n)
        else
          push!(topright_nodes, n)
        end
      end
    end

    if length(topleft_nodes) > 0
      if length(topleft_nodes) < length(r.nodes)
        subregion = Region(topleft_nodes)
        push!(r.subregions, subregion)
      else
        for n in topleft_nodes
          subregion = Region([n])
          push!(r.subregions, subregion)
        end
      end
    end

    if length(bottomleft_nodes) > 0
      if length(bottomleft_nodes) < length(r.nodes)
        subregion = Region(bottomleft_nodes)
        push!(r.subregions, subregion)
      else
        for n in bottomleft_nodes
          subregion = Region([n])
          push!(r.subregions, subregion)
        end
      end
    end

    if length(topright_nodes) > 0
      if length(topright_nodes) < length(r.nodes)
        subregion = Region(topright_nodes)
        push!(r.subregions, subregion)
      else
        for n in topright_nodes
          subregion = Region([n])
          push!(r.subregions, subregion)
        end
      end
    end

    if length(bottomright_nodes) > 0
      if length(bottomright_nodes) < length(r.nodes)
        subregion = Region(bottomright_nodes)
        push!(r.subregions, subregion)
      else
        for n in bottomright_nodes
          subregion = Region([n])
          push!(r.subregions, subregion)
        end
      end
    end

    for subregion in r.subregions
      build_sub_regions!(subregion)
    end
  end
end

function apply_force!(r::Region, n::Node, theta::Real; coefficient::Real = 0)
  if length(r.nodes) < 2
    linear_repulsion!(n, r.nodes[1]; coefficient = coefficient)
  else
    distance = sqrt((n.x - r.mass_center_x)^2 + (n.y - r.mass_center_y)^2)
    if distance * theta > r.size
      linear_repulsion_region!(n, r; coefficient = coefficient)
    else
      for subregion in r.subregions
        apply_force!(subregion, n, theta; coefficient = coefficient)
      end
    end
  end
end

function apply_force_on_nodes!(r::Region, nodes::Vector{Node}, theta::Real; coefficient::Real = 0)
  for n in nodes
    apply_force!(r, n, theta; coefficient = coefficient)
  end
end


Base.@kwdef struct AdjustSpeedAndApplyForcesResult
  speed::Real
  speed_efficiency::Real
end

function adjust_speed_and_apply_forces!(nodes::Vector{Node}, speed::Real, speed_efficiency::Real, jitter_tolerance::Real)::AdjustSpeedAndApplyForcesResult
  total_swinging = 0.0
  total_effective_traction = 0.0

  for n in nodes
    swinging = sqrt((n.old_dx - n.dx)^2 + (n.old_dy - n.dy)^2)
    total_swinging += n.mass * swinging
    total_effective_traction += 0.5 * n.mass * sqrt((n.old_dx + n.dx) * (n.old_dx + n.dx) + (n.old_dy + n.dy) * (n.old_dy + n.dy))
  end

  estimated_optinal_jitter_tolerance = 0.05 * sqrt(length(nodes))
  min_jt = sqrt(estimated_optinal_jitter_tolerance)
  max_jt = 10
  jt = jitter_tolerance * max(
    min_jt,
    min(max_jt, estimated_optinal_jitter_tolerance * total_effective_traction / (length(nodes)^2)),
  )
  min_speed_efficiency = 0.05

  if total_effective_traction > 0 && total_swinging / total_effective_traction > 2
    if speed_efficiency > min_speed_efficiency
      speed_efficiency *= 0.5
    end
    jt = max(jt, jitter_tolerance)
  end

  target_speed::Real = 0
  if total_swinging == 0
    target_speed = Inf
  else
    target_speed = jt * speed_efficiency * total_effective_traction / total_swinging
  end

  if total_swinging > jt * total_effective_traction
    if speed_efficiency > min_speed_efficiency
      speed_efficiency *= 0.7
    end
  elseif speed < 1000
    speed_efficiency *= 1.3
  end

  max_rise = 0.5
  speed = speed + min(target_speed - speed, max_rise * speed)

  for n in nodes
    swinging = n.mass * sqrt((n.old_dx - n.dx)^2 + (n.old_dy - n.dy)^2)
    factor = speed / (1.0 + sqrt(speed * swinging))
    n.x += n.dx * factor
    n.y += n.dy * factor
  end

  return AdjustSpeedAndApplyForcesResult(speed, speed_efficiency)
end