/*
  Copyright 2019 Balázs Ludmány

  This file is part of libcontours.

  libcontours is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  libcontours is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with libcontours.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef DISCOVER_GRAPH_HPP
#define DISCOVER_GRAPH_HPP 1

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/reverse_graph.hpp>
#include <boost/property_map/property_map.hpp>

namespace contours {
  template <typename Graph, typename VertexArea, typename EdgeArea, typename RootsMap>
static void
discover_from(const Graph &graph,
              const typename boost::graph_traits<Graph>::vertex_descriptor root,
	      const VertexArea &vertex_area,
              EdgeArea &edge_area, RootsMap &roots) {
  using namespace boost;
  using Real = typename property_traits<EdgeArea>::value_type;
  using Set = typename property_traits<RootsMap>::value_type;

  if (out_degree(root, graph) != 1)
    return;

  Real sum_area = get(vertex_area, root);
  Set sum_roots;
  for (const auto edge : make_iterator_range(in_edges(root, graph))) {
    if (get(edge_area, edge) == 0.0)
      return;
    sum_area += get(edge_area, edge);
    sum_roots.merge(get(roots, edge));
  }
  if (sum_roots.empty())
    sum_roots.insert(root);

  const auto next = *out_edges(root, graph).first;
  put(edge_area, next, sum_area);
  put(roots, next, sum_roots);

  discover_from(graph, target(next, graph), vertex_area, edge_area, roots);
}

  template <typename Graph, typename VertexArea, typename EdgeArea, typename RootsMap,
          typename VertexIterator>
void discover_graph(const Graph &graph, const VertexArea &vertex_area, EdgeArea &edge_area, RootsMap &roots, VertexIterator leaves) {
  using namespace boost;

  for (const auto vertex : make_iterator_range(vertices(graph))) {
    if (in_degree(vertex, graph) == 0 && out_degree(vertex, graph) != 0) {
      discover_from(graph, vertex, vertex_area, edge_area, roots);
      *(leaves++) = vertex;
    }
  }
}
}; // namespace shape

#endif // DISCOVER_GRAPH_HPP
